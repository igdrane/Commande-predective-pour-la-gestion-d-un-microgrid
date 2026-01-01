%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbre_jours=7;
Te=10;
Hp=24*60/Te;
t=0:Te:Nbre_jours*24*60-Te;
N_fin=length(t);
T_fin=Te*N_fin;
%le vecteur des tarifs en euro/KWh:
T_HC=0.1228;
T_HP=0.1579;
T_EDF_jour = [ ...
    T_HC*ones(6*60/Te,1); ...
    T_HP*ones((8*60-6*60)/Te,1); ...
    T_HC*ones((12*60-8*60)/Te,1); ...
    T_HP*ones((14*60-12*60)/Te,1); ...
    T_HC*ones((16*60-14*60)/Te,1); ...
    T_HP*ones((22*60-16*60)/Te,1); ...
    T_HC*ones((24*60-22*60)/Te,1) ...
];
T_EDF=[];
for i=1:Nbre_jours
     T_EDF=[T_EDF; T_EDF_jour];
end
load('data_exemple.mat');

%La production photovoltaique:
S=40;
r=0.2;
Cp=0.1;
Pp = S * r * GHI_real * (1 - Cp);
Pp_pred = S * r * GHI_cs * (1 - Cp);
%Les parametres de la batterie:
Cb=100;
Eb0=0.1*Cb;
Eb(1)=Eb0;
Pb(1)=0;

%Calcule de Eb(2) et Pbr(1) et Prr:
[Eb(2), Pbr(1)]=systeme(Pb(1),Eb(1),Cb,Te);
Prr(1)=Pc(1)+Pbr(1)-Pp(1);

%% 47. Graphiques de comparaison des prédictions
% Graphique 1: GHI réel vs prédit
subplot(3,1,1);
plot(temps_data, GHI_real, 'b-', 'LineWidth', 2, 'DisplayName', 'GHI réel');
hold on;
plot(temps_data, GHI_cs, 'r--', 'LineWidth', 1.5, 'DisplayName', 'GHI prédit (ciel clair)');
xlabel('Temps');
ylabel('GHI (kW/m²)');
title('Comparaison du GHI réel et prédit');
legend('Location', 'best');
grid on;

% Graphique 2: Production PV réelle vs prédite
subplot(3,1,2);

plot(temps_data, Pp, 'b-', 'LineWidth', 2, 'DisplayName', 'Production réelle');
hold on;
plot(temps_data, Pp_pred, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Production prédite');
xlabel('Temps');
ylabel('Puissance (kW)');
title('Comparaison de la production PV réelle et prédite');
legend('Location', 'best');
grid on;

% Graphique 3: Consommation réelle vs prédite
subplot(3,1,3);
plot(temps_data, Pc, 'b-', 'LineWidth', 2, 'DisplayName', 'Consommation réelle');
hold on;
plot(temps_data, Pc_pred, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Consommation prédite');
xlabel('Temps');
ylabel('Puissance (kW)');
title('Comparaison de la consommation réelle et prédite');
legend('Location', 'best');
grid on;

%% 48. Calcul des fits et analyse
fprintf('\n--- 48. Calcul des fits des prédictions ---\n');

% Calcul des fits
fit_production = calcul_fit(Pp, Pp_pred);
fit_consommation = calcul_fit(Pc, Pc_pred);
fit_GHI=calcul_fit(GHI_real, GHI_cs);


fprintf('Fit de GHI: %.1f%%\n', fit_GHI);
fprintf('Fit de la production PV: %.1f%%\n', fit_production);
fprintf('Fit de la consommation: %.1f%%\n', fit_consommation);


%% Boucle de simulation et d:optimisation en temps rée: Etude des cas:
%Matrices pour l'enregstriment des differents parametres suivants les cas:
M_Eb=zeros(1009, 3);
M_Pbr=zeros(1008, 3);
M_Pr=zeros(1008, 3);
M_Prr=zeros(1008, 3);

%les cas de simulations:
Cas_data=[ "Pc", "Pp_pred"; "Pc_pred", "Pp"; "Pc_pred", "Pp_pred"];

for i=1:3
    switch Cas_data(i,:)
    case ["Pc", "Pp_pred"]
        Pc_utilise=Pc;
        Pp_utilise=Pp_pred;

        % Pr(1)=Pc(1)+Pb(1)-Pp_pred(1);
        % [Eb(2), Pbr(1)]=systeme(Pb(1),Eb(1),Cb,Te);
        % Prr(1)=Pc(1)+Pbr(1)-Pp_pred(1);
    case ["Pc_pred", "Pp"]
        Pc_utilise=Pc_pred;
        Pp_utilise=Pp;

        % Pr(1)=Pc_pred(1)+Pb(1)-Pp(1);
        % [Eb(2), Pbr(1)]=systeme(Pb(1),Eb(1),Cb,Te);
        % Prr(1)=Pc_pred(1)+Pbr(1)-Pp(1);
    case ["Pc_pred", "Pp_pred"]
        Pc_utilise=Pc_pred;
        Pp_utilise=Pp_pred;

        % Pr(1)=Pc_pred(1)+Pb(1)-Pp_pred(1);
        % [Eb(2), Pbr(1)]=systeme(Pb(1),Eb(1),Cb,Te);
        % Prr(1)=Pc_pred(1)+Pbr(1)-Pp_pred(1);
    end
    for k=1:N_fin-1
        disp(['cas: ', num2str(i), ' itération: ', num2str(k), '/', num2str(N_fin-1)]);
        Hpr=min(Hp, N_fin-k);
        T_EDF_H=T_EDF(k+1:k+Hpr);
        Pr0=ones(Hpr,1);
        Prmin=-9*ones(Hpr,1);
        Prmax=9*ones(Hpr,1);
        Pc_H=Pc_utilise(k+1:k+Hpr);
        Pp_H=Pp_utilise(k+1:k+Hpr);
        options=optimoptions("fmincon","Algorithm","sqp","Display","off","MaxFunctionEvaluations",10000,"MaxIterations",1000);
        fob=@(Pr) funobj(Pr, Eb(k+1), T_EDF_H, Pc_H, Pp_H, Cb, Prmax, Prmin, Te);
        [Pr_opt, fob_val, exitflag]=fmincon(fob,Pr0,[],[],[],[],Prmin,Prmax,[],options);
        Pr(k+1)=Pr_opt(1);
        Pb(k+1)=Pp(k+1)+Pr(k+1)-Pc(k+1);
        [Eb(k+2), Pbr(k+1)]=systeme(Pb(k+1), Eb(k+1), Cb, Te);
        Prr(k+1)=Pc(k+1)+Pbr(k+1)-Pp(k+1);
    end

    M_Eb(:, i)=Eb';
   
    M_Pbr(:,i)=Pbr';
   
    M_Pr(:,i)=Pr';
   
    M_Prr(:,i)=Prr';
    
    Eb(1)=Eb0;
    Pb(1)=0;
    Pr(1)=Pc(1)+Pb(1)-Pp(1);
    [Eb(2), Pbr(1)]=systeme(Pb(1),Eb(1),Cb,Te);
    Prr(1)=Pc(1)+Pbr(1)-Pp(1);
end

save data_resultats_cas M_Eb M_Pbr M_Pr M_Prr
%le vecteut cout Cr de l'electricité echangé avec le reseau:
Cr_cas1=(M_Prr(:,1)).*T_EDF*Te/60;
Cr_cas2=(M_Prr(:,2)).*T_EDF*Te/60;
Cr_cas3=(M_Prr(:,3)).*T_EDF*Te/60;

%le cout Cr de l'electricité echangé avec le reseau:
Cr_tot_cas1=sum(Cr_cas1);
Cr_tot_cas2=sum(Cr_cas2);
Cr_tot_cas3=sum(Cr_cas3);

fprintf("le cout total cas 1 est: %f ", Cr_tot_cas1)
sprintf("le cout total cas 2 est: %f ", Cr_tot_cas2)
sprintf("le cout total cas 3 est: %f ", Cr_tot_cas3)




%% l'affichage des resultats cas 1: "Pc", "Pp_pred"
figure
subplot(3,1,1)
tb=[t T_fin+Te];
xlim([0, tb(end)])
title('Evolution Energie de la batterie')
grid
hold on
Emin=10;
Emax=100;
plot(tb,M_Eb(:,1))
xlabel('t en minutes')
ylabel('Energie en KWh')
plot(tb,Emin*ones(length(tb),1))
plot(tb,Emax*ones(length(tb),1))
legend('Eb','Emin','Emax')
hold off

subplot(3,1,2)
xlim([0, t(end)])
title('Evolution de puissance Pc, Pp, Pbr et Prr')
grid
hold on
plot(t, Pc);
plot(t, Pp_pred);
plot(t, M_Pbr(:,1));
plot(t, M_Prr(:,1));
plot(t,-9*ones(length(t),1))
plot(t,9*ones(length(t),1))
legend('Pc','Pp-pred','Pbr','Prr','Pmin','Pmax')
xlabel('t en minutes')
ylabel('Puissance en KW')
hold off

subplot(3,1,3)
xlim([0, t(end)])
title('Le cout total et les tarifs')
grid
hold on
yyaxis left;
plot(t, Cr_cas1)
ylabel('Cout en euro')
yyaxis right;
plot(t, T_EDF)
xlabel('t en minutes')
ylabel('Tarif en euro/KWh')
hold off
sgtitle("l'affichage des resultats cas 1: Pc, Pp-pred")


%% l'affichage des resultats cas 2: "Pc_pred", "Pp"
figure
subplot(3,1,1)
tb=[t T_fin+Te];
xlim([0, tb(end)])
title('Evolution Energie de la batterie')
grid
hold on
Emin=10;
Emax=100;
plot(tb,M_Eb(:,2))
xlabel('t en minutes')
ylabel('Energie en KWh')
plot(tb,Emin*ones(length(tb),1))
plot(tb,Emax*ones(length(tb),1))
legend('Eb','Emin','Emax')
hold off

subplot(3,1,2)
xlim([0, t(end)])
title('Evolution de puissance Pc, Pp, Pbr et Prr')
grid
hold on
plot(t, Pc_pred);
plot(t, Pp);
plot(t, M_Pbr(:,2));
plot(t, M_Prr(:,2));
plot(t,-9*ones(length(t),1))
plot(t,9*ones(length(t),1))
legend('Pc-pred','Pp','Pbr','Prr','Pmin','Pmax')
xlabel('t en minutes')
ylabel('Puissance en KW')
hold off

subplot(3,1,3)
xlim([0, t(end)])
title('Le cout total et les tarifs')
grid
hold on
yyaxis left;
plot(t, Cr_cas2)
ylabel('Cout en euro')
yyaxis right;
plot(t, T_EDF)
xlabel('t en minutes')
ylabel('Tarif en euro/KWh')
hold off
sgtitle("l'affichage des resultats cas 2: Pc-pred, Pp")




%% l'affichage des resultats cas 3: "Pc_pred", "Pp_pred"
figure
subplot(3,1,1)
tb=[t T_fin+Te];
xlim([0, tb(end)])
title('Evolution Energie de la batterie')
grid
hold on
Emin=10;
Emax=100;
plot(tb,M_Eb(:,3))
xlabel('t en minutes')
ylabel('Energie en KWh')
plot(tb,Emin*ones(length(tb),1))
plot(tb,Emax*ones(length(tb),1))
legend('Eb','Emin','Emax')
hold off

subplot(3,1,2)
xlim([0, t(end)])
title('Evolution de puissance Pc, Pp, Pbr et Prr')
grid
hold on
plot(t, Pc_pred);
plot(t, Pp_pred);
plot(t, M_Pbr(:,3));
plot(t, M_Prr(:,3));
plot(t,-9*ones(length(t),1))
plot(t,9*ones(length(t),1))
legend('Pc-pred','Pp-pred','Pbr','Prr','Pmin','Pmax')
xlabel('t en minutes')
ylabel('Puissance en KW')
hold off

subplot(3,1,3)
xlim([0, t(end)])
title('Le cout total et les tarifs')
grid
hold on
yyaxis left;
plot(t, Cr_cas3)
ylabel('Cout en euro')
yyaxis right;
plot(t, T_EDF)
xlabel('t en minutes')
ylabel('Tarif en euro/KWh')
hold off
sgtitle("l'affichage des resultats cas 3: Pc-pred, Pp-pred")










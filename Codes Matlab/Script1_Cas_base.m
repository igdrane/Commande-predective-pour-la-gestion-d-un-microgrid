clear
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbre_jours=7;
Te=10;
t=0:Te:Nbre_jours*24*60-Te;
N_fin=length(t);
T_fin=N_fin*Te;
%Horizone de prediction:
Hp=24*60/Te;
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

%Données
load('data_exemple.mat');

%La production photovoltaique:
S=40;
r=0.2;
Cp=0.1;
Pp=S*r*GHI_real*(1-Cp);

%Les parametres de la batterie:
Cb=100;
Eb0=0.1*Cb;
Eb(1)=Eb0;
Pb(1)=0;


Pr(1)=Pc(1)+Pb(1)-Pp(1);
%Calcule de Eb(2) et Pbr(1):
[Eb(2), Pbr(1)]=systeme(Pb(1),Eb(1),Cb,Te);

Prr(1)=Pc(1)+Pbr(1)-Pp(1);
fprintf('la valeur de Eb(2) est: %.3f KW\n', Eb(2))
fprintf('la valeur de Pbr(1) est: %.3f KW\n', Pbr(1))
fprintf('la valeur de Prr(1) est: %.3f KW\n', Prr(1))


options=optimoptions("fmincon","Algorithm","sqp","Display","off","MaxFunctionEvaluations",10000,"MaxIterations",1000);
%Partie 2.3.2 Boucle de simulation et d’optimisation en temps réel 
for k=1:N_fin-1
    disp(['itération: ', num2str(k), '/', num2str(N_fin-1)]);
    Hpr=min(Hp, N_fin-k);
    T_EDF_H=T_EDF(k+1:k+Hpr);
    Pc_H=Pc(k+1:k+Hpr);
    Pp_H=Pp(k+1:k+Hpr);
    Pr0=ones(Hpr,1);
    Prmin=-9*ones(Hpr,1);
    Prmax=9*ones(Hpr,1);

    fob=@(Pr) funobj(Pr, Eb(k+1), T_EDF_H, Pc_H, Pp_H, Cb, Prmax, Prmin, Te);
    [Pr_opt, fob_val, exitflag]=fmincon(fob,Pr0,[],[],[],[],Prmin,Prmax,[],options);
    Pr(k+1)=Pr_opt(1);
    Pb(k+1)=Pp(k+1)+Pr(k+1)-Pc(k+1);
    [Eb(k+2), Pbr(k+1)]=systeme(Pb(k+1), Eb(k+1), Cb, Te);
    Prr(k+1)=Pc(k+1)+Pbr(k+1)-Pp(k+1);

    disp(['fob = ', num2str(fob_val)]);
    disp(['exitflag = ', num2str(exitflag)]); 
end 

%le vecteut cout Cr de l'electricité echangé avec le reseau:
Cr=(Prr)'.*T_EDF*Te/60;

%le cout Cr de l'electricité echangé avec le reseau:
Cr_tot=sum(Cr);
sprintf("le cout total est: %f ", Cr_tot)

%l'affichage des resultats:
subplot(3,1,1)
tb=[t T_fin+Te];
xlim([0, tb(end)])
title('Evolution Energie de la batterie')
grid
hold on
Emin=10;
Emax=100;
plot(tb,Eb)
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
plot(t, Pp);
plot(t, Pbr);
plot(t, Prr);
plot(tb,-9*ones(length(tb),1))
plot(tb,9*ones(length(tb),1))
legend('Pc','Pp','Pbr','Prr','Pmin','Pmax')
xlabel('t en minutes')
ylabel('Puissance en KW')
hold off

subplot(3,1,3)
xlim([0, t(end)])
title('Le cout total et les tarifs')
grid
hold on
yyaxis left;
plot(t, Cr)
ylabel('Cout en euro')
yyaxis right;
plot(t, T_EDF)
xlabel('t en minutes')
ylabel('Tarif en euro/KWh')
hold off














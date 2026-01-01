clear
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMÈTRES
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbre_jours = 7;
Te = 10;
Hp = 24*60/Te;

t = 0:Te:Nbre_jours*24*60-Te;
N_fin = length(t);
T_fin = N_fin*Te;

% Tarifs
T_HC = 0.1228;
T_HP = 0.1579;
T_EDF_jour = [ ...
    T_HC*ones(6*60/Te,1); ...
    T_HP*ones((8*60-6*60)/Te,1); ...
    T_HC*ones((12*60-8*60)/Te,1); ...
    T_HP*ones((14*60-12*60)/Te,1); ...
    T_HC*ones((16*60-14*60)/Te,1); ...
    T_HP*ones((22*60-16*60)/Te,1); ...
    T_HC*ones((24*60-22*60)/Te,1) ...
];

T_EDF = [];
for i=1:Nbre_jours
    T_EDF = [T_EDF; T_EDF_jour];
end

% Données
load('data_exemple.mat');   % Pc, GHI_real, ...

% PV
S=40; r=0.2; Cp=0.1;
Pp_pred = S * r * GHI_cs * (1 - Cp);
Pp = S*r*GHI_real*(1-Cp);

% Batterie
Cb=100;
Eb0=0.1*Cb;

% Options fmincon
options = optimoptions("fmincon", ...
    "Algorithm","sqp","Display","off", ...
    "MaxFunctionEvaluations",10000, ...
    "MaxIterations",1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 INITIALISATIONS À TESTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%
cas_Pr0 = ["Pr0=1", "Pr0=0", "Pr0=Pc_H-Pp_H"];

% Stockage résultats
M_Eb  = zeros(N_fin+1, 3);
M_Pbr = zeros(N_fin,   3);
M_Pr  = zeros(N_fin,   3);
M_Prr = zeros(N_fin,   3);
Cr_tot = zeros(3,1);
Time_tot = zeros(3,1);

for i = 1:3
    % Réinitialisation variables
    Eb  = zeros(N_fin+1,1);
    Pb  = zeros(N_fin,1);
    Pr  = zeros(N_fin,1);
    Pbr = zeros(N_fin,1);
    Prr = zeros(N_fin,1);

    Eb(1)=Eb0;
    Pb(1)=0;
    Pr(1)=Pc(1)+Pb(1)-Pp(1);

    [Eb(2), Pbr(1)] = systeme(Pb(1), Eb(1), Cb, Te);
    Prr(1)=Pc(1)+Pbr(1)-Pp(1);

    t_run = tic;

    for k=1:N_fin-1
        disp(['itération: ', num2str(k), '/', num2str(N_fin-1)]);

        Hpr = min(Hp, N_fin-k);
        T_EDF_H = T_EDF(k+1:k+Hpr);
        Pc_H = Pc_pred(k+1:k+Hpr);
        Pp_H = Pp_pred(k+1:k+Hpr);

        Prmin = -9*ones(Hpr,1);
        Prmax = 9*ones(Hpr,1);

        % --------- choix de l'initialisation Pr0 ----------
        switch cas_Pr0(i)
            case "Pr0=1"
                Pr0 = ones(Hpr,1);
            case "Pr0=0"
                Pr0 = zeros(Hpr,1);
            case "Pr0=Pc_H-Pp_H"
                Pr0 = Pc_H - Pp_H;
        end
        % ---------------------------------------------------

        fob = @(Pr) funobj(Pr, Eb(k+1), T_EDF_H, Pc_H, Pp_H, Cb, Prmax, Prmin, Te);
        [Pr_opt, fob_val, exitflag] = fmincon(fob, Pr0, [],[],[],[], Prmin, Prmax, [], options);

        Pr(k+1)=Pr_opt(1);

        % Simulation réelle
        Pb(k+1)=Pp(k+1)+Pr(k+1)-Pc(k+1);
        [Eb(k+2), Pbr(k+1)] = systeme(Pb(k+1), Eb(k+1), Cb, Te);
        Prr(k+1)=Pc(k+1)+Pbr(k+1)-Pp(k+1);
    end

    Time_tot(i) = toc(t_run);

    % Coût total
    Cr = (Prr).*T_EDF*Te/60;
    Cr_tot(i) = sum(Cr);

    % Stockage
    M_Eb(:,i)  = Eb;
    M_Pbr(:,i) = Pbr;
    M_Pr(:,i)  = Pr;
    M_Prr(:,i) = Prr;

    fprintf(">> %s : Cout total = %.4f € | Temps total = %.3f s\n", cas_Pr0(i), Cr_tot(i), Time_tot(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARAISONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(3,1,1); grid on; hold on;
tb = [t, t(end)+Te];
plot(tb, M_Eb(:,1));
plot(tb, M_Eb(:,2));
plot(tb, M_Eb(:,3));
xlabel('t (min)'); ylabel('Eb (kWh)');
title('Eb : comparaison des 3 initialisations');
legend(cas_Pr0(:), 'Location','best');

subplot(3,1,2); grid on; hold on;
plot(t, M_Prr(:,1));
plot(t, M_Prr(:,2));
plot(t, M_Prr(:,3));
plot(t, -9*ones(size(t)));
plot(t, 9*ones(size(t)));
xlabel('t (min)'); ylabel('Prr (kW)');
title('Prr : comparaison');
legend([cas_Pr0, {'Pmin','Pmax'}], 'Location','best');

subplot(3,1,3); grid on;
bar(Cr_tot);
set(gca,'XTickLabel', cas_Pr0);
ylabel('Coût total (€)');
title('Coût total selon Pr0');

figure;
bar(Time_tot);
set(gca,'XTickLabel', cas_Pr0);
ylabel('Temps total (s)');
title('Temps de calcul total selon Pr0');

% save resultats_3init_Pr0 M_Eb M_Pbr M_Pr M_Prr Cr_tot Time_tot cas_Pr0

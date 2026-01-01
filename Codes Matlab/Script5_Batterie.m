clear; clc; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMÈTRES (comme ton script)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nbre_jours = 7;
Te = 10;
Hp = 24*60/Te;

t = 0:Te:Nbre_jours*24*60-Te;
N_fin = length(t);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DONNÉES + PV
%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('data_exemple.mat');  % doit contenir Pc, GHI_real

S=40; r=0.2; Cp=0.1;
Pp = S*r*GHI_real*(1-Cp);
Pp_pred=S*r*GHI_cs*(1-Cp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVEUR (SQP comme ton script)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Prmin_val = -9;
Prmax_val =  9;

opt_pso=optimoptions("particleswarm","Display","off", ...
    "SwarmSize",200,"MaxIterations",1000,"FunctionTolerance",1e-6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ÉTUDE CAPACITÉ: 0 à 100 kWh (pas de 10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cb_list = (0:10:100)';
nCb = length(Cb_list);

Cout_7j   = zeros(nCb,1);   % coût sur 7 jours
Cost_10a  = zeros(nCb,1);   % coût extrapolé sur 10 ans
Gain_7j   = zeros(nCb,1);   % gain brut sur 7 jours (vs sans batterie)
Gain_10a  = zeros(nCb,1);   % gain brut sur 10 ans
Net_10a_2019 = zeros(nCb,1);% gain net 10 ans - prix batterie (200 €/kWh)
Net_10a_2030 = zeros(nCb,1);% gain net 10 ans - prix batterie (100 €/kWh)

% ---------------------------
% Calcul coût sans batterie (Cb=0) -> référence
% ---------------------------
Prr0 = Pc - Pp;                       % Puissance reseau sans batterie
Cr0  = (Prr0).*T_EDF*Te/60;
Cout_7d_0 = sum(Cr0);

for i = 1:nCb
    Cb = Cb_list(i);

    if Cb == 0
        % Sans batterie : déjà calculé
        Cout_7j(i) = Cout_7d_0;

    else
        % ---------------------------
        % Avec batterie
        % ---------------------------
        Eb0 = 0.1*Cb;

        Eb  = zeros(1, N_fin+1);
        Pb  = zeros(1, N_fin);
        Pr  = zeros(1, N_fin);
        Pbr = zeros(1, N_fin);
        Prr = zeros(1, N_fin);

        Eb(1) = Eb0;
        Pb(1) = 0;

        Pr(1) = Pc(1) + Pb(1) - Pp(1);
        [Eb(2), Pbr(1)] = systeme(Pb(1), Eb(1), Cb, Te);
        Prr(1) = Pc(1) + Pbr(1) - Pp(1);

        for k = 1:N_fin-1
            disp(['Cas Cb= ',num2str(Cb) ,'itération: ', num2str(k), '/', num2str(N_fin-1)]);

            Hpr = min(Hp, N_fin-k);

            T_EDF_h = T_EDF(k+1:k+Hpr);
            Pc_h    = Pc_pred(k+1:k+Hpr);
            Pp_h    = Pp_pred(k+1:k+Hpr);

            Prmin = -9*ones(Hpr,1);
            Prmax = 9*ones(Hpr,1);

            Pr0 = ones(Hpr,1);

            fob = @(Pr) funobj(Pr, Eb(k+1), T_EDF_h, Pc_h, Pp_h, Cb, Prmax, Prmin, Te);
            [Pr_opt, fob_val, exitflag]=particleswarm(fob,Hpr,Prmin,Prmax,opt_pso);


            Pr(k+1) = Pr_opt(1);

            Pb(k+1) = Pp(k+1) + Pr(k+1) - Pc(k+1);
            [Eb(k+2), Pbr(k+1)] = systeme(Pb(k+1), Eb(k+1), Cb, Te);
            Prr(k+1) = Pc(k+1) + Pbr(k+1) - Pp(k+1);
        end

        Cr = (Prr)'.*T_EDF*Te/60;
        Cout_7j(i) = sum(Cr);
    end

    % Gains / extrapolations
    Gain_7j(i)  = Cout_7d_0 - Cout_7j(i);          % positif = économie
    Gain_10a(i) = Gain_7j(i) * 52 * 10;
    Cost_10a(i) = Cout_7j(i) * 52 * 10;

    Net_10a_2019(i) = Gain_10a(i) - 200*Cb;        % 200 €/kWh
    Net_10a_2030(i) = Gain_10a(i) - 100*Cb;        % 100 €/kWh
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 53-54) Coûts (Cb=0 et Cb=100)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, i0]   = ismember(0, Cb_list);
[~, i100] = ismember(100, Cb_list);

fprintf("\n53) Coût total sans batterie (Cb=0 kWh) sur 7 jours : %.4f €\n", Cout_7j(i0));
fprintf("54) Coût total avec Cb=100 kWh sur 7 jours : %.4f €\n", Cout_7j(i100));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 55) Capacité maximale utile
% Définition utilisée : plus petite capacité atteignant 99%% du gain à 100 kWh
%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain_ref = Gain_7j(i100);
seuil = 0.99 * gain_ref;

idx_util = find(Gain_7j >= seuil, 1, 'first');
Cb_utile = Cb_list(idx_util);

fprintf("\n55) Capacité maximale utile (>=99%% du gain à 100 kWh) : %d kWh\n", Cb_utile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 56-58) Gains bruts (hebdo, annuel, 10 ans)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
gain7_utile  = Gain_7j(idx_util);
gain52_utile = gain7_utile * 52;
gain10y_utile = gain7_utile * 52 * 10;

fprintf("\n56) Gain brut sur 7 jours avec Cb_utile : %.4f €\n", gain7_utile);
fprintf("57) Gain brut extrapolé sur 52 semaines : %.4f € / an\n", gain52_utile);
fprintf("58) Gain brut extrapolé sur 10 ans : %.4f €\n", gain10y_utile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 59-60) Prix 2019 (200 €/kWh) + gain/perte net sur 10 ans
%%%%%%%%%%%%%%%%%%%%%%%%%%%
prix2019 = 200 * Cb_utile;
net10y_2019 = gain10y_utile - prix2019;

fprintf("\n59) Prix batterie (2019, 200 €/kWh) avec Cb_utile : %.2f €\n", prix2019);
fprintf("60) Gain/Pertes net sur 10 ans (2019) : %.4f €\n", net10y_2019);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 61-62) Prix 2030 (100 €/kWh) + gain/perte net sur 10 ans
%%%%%%%%%%%%%%%%%%%%%%%%%%%
prix2030 = 100 * Cb_utile;
net10y_2030 = gain10y_utile - prix2030;

fprintf("\n61) Prix batterie (2030, 100 €/kWh) avec Cb_utile : %.2f €\n", prix2030);
fprintf("62) Gain/Pertes net sur 10 ans (2030) : %.4f €\n", net10y_2030);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 63) Graph coût total extrapolé sur 10 ans vs capacité (pas 10 kWh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(Cb_list, Cost_10a, '-o');
grid on;
xlabel('Capacité batterie Cb (kWh)');
ylabel('Coût total extrapolé sur 10 ans (€)');
title('63) Coût total (10 ans) en fonction de la capacité batterie');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 64) Graph gain brut extrapolé sur 10 ans vs capacité
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(Cb_list, Gain_10a, '-o');
grid on;
xlabel('Capacité batterie Cb (kWh)');
ylabel('Gain brut extrapolé sur 10 ans (€)');
title('64) Gain brut (10 ans) en fonction de la capacité batterie');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 65) Graph gain/perte net sur 10 ans incluant prix d’achat
% (On trace les deux scénarios 2019 et 2030)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(Cb_list, Net_10a_2019, '-o'); hold on;
plot(Cb_list, Net_10a_2030, '-s');
yline(0,'--');
grid on;
xlabel('Capacité batterie Cb (kWh)');
ylabel('Gain/Pertes net sur 10 ans (€)');
title('65) Gain/Pertes net (10 ans) incluant prix batterie');
legend('2019 (200 €/kWh)','2030 (100 €/kWh)','Location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 66) Recommandation (simple et basée sur le gain net max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
[netmax2019, idxm2019] = max(Net_10a_2019);
[netmax2030, idxm2030] = max(Net_10a_2030);

fprintf("\n66) Recommandation:\n");
fprintf(" - Optimum net 2019 : Cb = %d kWh, gain net 10 ans = %.2f €\n", Cb_list(idxm2019), netmax2019);
fprintf(" - Optimum net 2030 : Cb = %d kWh, gain net 10 ans = %.2f €\n", Cb_list(idxm2030), netmax2030);

if netmax2019 <= 0
    fprintf("   -> En 2019 (200 €/kWh), la batterie n’est pas rentable sur 10 ans (selon ces données).\n");
else
    fprintf("   -> En 2019 (200 €/kWh), une capacité proche de %d kWh est recommandée.\n", Cb_list(idxm2019));
end

if netmax2030 <= 0
    fprintf("   -> En 2030 (100 €/kWh), la batterie n’est pas rentable sur 10 ans (selon ces données).\n");
else
    fprintf("   -> En 2030 (100 €/kWh), une capacité proche de %d kWh est recommandée.\n", Cb_list(idxm2030));
end

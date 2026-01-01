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
Pp_pred=S*r*GHI_cs*(1-Cp);

%Les parametres de la batterie:
Cb=100;
Eb0=0.1*Cb;

%Bornes
Prmin_val = -9;
Prmax_val =  9;

%Options SQP
opt_sqp=optimoptions("fmincon","Algorithm","sqp","Display","off", ...
    "MaxFunctionEvaluations",10000,"MaxIterations",1000);

%Options PSO
opt_pso=optimoptions("particleswarm","Display","off", ...
    "SwarmSize",100,"MaxIterations",1000,"FunctionTolerance",1e-6);

%On compare les 2 algorithmes
algo_list = ["SQP","PSO"];
Cr_tot_all = zeros(2,1);
Time_all   = zeros(2,1);

for a=1:2
    algo = algo_list(a);

    % ===========================
    % INITIALISATION (comme base)
    % ===========================
    Eb=zeros(1,N_fin+1);
    Pb=zeros(1,N_fin);
    Pr=zeros(1,N_fin);
    Pbr=zeros(1,N_fin);
    Prr=zeros(1,N_fin);

    Eb(1)=Eb0;
    Pb(1)=0;

    Pr(1)=Pc(1)+Pb(1)-Pp(1);
    [Eb(2), Pbr(1)]=systeme(Pb(1),Eb(1),Cb,Te);
    Prr(1)=Pc(1)+Pbr(1)-Pp(1);

    fprintf("\n========================\n");
    fprintf("Algorithme : %s\n", algo);
    fprintf("========================\n");

    t_run = tic;

    % ===========================
    % BOUCLE MPC (comme base)
    % ===========================
    for k=1:N_fin-1

        disp(['itération: ', num2str(k), '/', num2str(N_fin-1)]);
        Hpr=min(Hp, N_fin-k);
        T_EDF_pred=T_EDF(k+1:k+Hpr);
        Pc_H=Pc_pred(k+1:k+Hpr);
        Pp_H=Pp_pred(k+1:k+Hpr);

        Prmin=Prmin_val*ones(Hpr,1);
        Prmax=Prmax_val*ones(Hpr,1);

        % point initial (tu peux garder ones comme dans ton base)
        Pr0=ones(Hpr,1);

        fob=@(Pr_var) funobj(Pr_var, Eb(k+1), T_EDF_pred, Pc_H, Pp_H, Cb, Prmax, Prmin, Te);

        switch algo
            case "SQP"
                [Pr_opt, fob_val, exitflag]=fmincon(fob,Pr0,[],[],[],[],Prmin,Prmax,[],opt_sqp);

            case "PSO"
                % particleswarm : nvars = Hpr
                [Pr_opt, fob_val, exitflag]=particleswarm(fob,Hpr,Prmin,Prmax,opt_pso);
        end

        Pr(k+1)=Pr_opt(1);
        Pb(k+1)=Pp(k+1)+Pr(k+1)-Pc(k+1);
        [Eb(k+2), Pbr(k+1)]=systeme(Pb(k+1), Eb(k+1), Cb, Te);
        Prr(k+1)=Pc(k+1)+Pbr(k+1)-Pp(k+1);
    end

    Time_all(a)=toc(t_run);

    % Cout total (comme base)
    Cr=(Prr)'.*T_EDF*Te/60;
    Cr_tot=sum(Cr);

    Cr_tot_all(a)=Cr_tot;

    fprintf(">> Cout total = %.4f\n", Cr_tot);
    fprintf(">> Temps total = %.3f s\n", Time_all(a));
end

fprintf("\n===== COMPARAISON FINALE =====\n");
fprintf("SQP : Cout = %.4f | Temps = %.3f s\n", Cr_tot_all(1), Time_all(1));
fprintf("PSO : Cout = %.4f | Temps = %.3f s\n", Cr_tot_all(2), Time_all(2));

% Choix meilleur (cout min)
[~, idx_best]=min(Cr_tot_all);
fprintf(">> Meilleur (cout min) : %s\n", algo_list(idx_best));

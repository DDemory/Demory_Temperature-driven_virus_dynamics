%% A thermal trade-off between viral production and degradation drives phytoplankton-virus population dynamics
% === Code for: ===
% SI Tab 6 -- Calcul of the cardinal parameters
% David Demory - Jan 2021

%% setup

T = [0:0.0001:33.5];
names = {'Mic-B/MicV-B','Mic-A/MicV-A','Mic-C/MicV-C'};
TTmin = []; TTopt = []; TTmax = []; Mopt = []; TToptR0 = []; TTR01 = [];
R0opt = [];

%% Calcul of the hyperparameters
for s = 1:3
    
    if s == 1;
        load('pmin_829_Arrhenius_v2.mat')
        [mu,K,phi,lambda,epsilon,beta,sigma,delta,omega,psi] = Tdriven_fct(T,pmin);
        [Topt,Tmax,muopt,Tmin] = Calculate_CT(pmin);
    elseif s == 2;
        %load('pmin_451_v3.mat')
        load('pmin_451_Arrhenius_v2.mat')
        [mu,K,phi,lambda,epsilon,beta,sigma,delta,omega,psi] = Tdriven_fct(T,pmin);
        [Topt,Tmax,muopt,Tmin] = Calculate_CT(pmin);
    elseif s == 3
        load('pmin_834_Arrhenius_v2.mat')
        [mu,K,phi,lambda,epsilon,beta,sigma,delta,omega,psi] = Tdriven_fct(T,pmin);
        [Topt,Tmax,muopt,Tmin] = Calculate_CT(pmin);
    end
    
    %% Calcul R0
    R0 = ((1-epsilon).*beta.*lambda.*phi.*K)./((lambda+psi).*(phi*K+sigma));
    ToptR0 = T(find(R0 == max(R0)));
    tempT = T(find(round(R0,2)==1));TR01 = tempT(1);

    %% Store Cardinal Parameters
    TTmin   = [TTmin,Tmin];
    TTopt   = [TTopt,Topt];
    TTmax   = [TTmax,Tmax];
    Mopt    = [Mopt,muopt];
    TToptR0 = [TToptR0,ToptR0];
    TTR01   = [TTR01,TR01];
    R0opt = [R0opt,max(R0)];
 
end

Cardinal_Parameters = {'Topt_mu','Tmax_mu','Tmin_mu','\mu_opt','Topt_R0','T_R01','R0_opt'}';
MicB_MicVB = [TTopt(1),TTmax(1),TTmin(1),Mopt(1),TToptR0(1),TTR01(1),R0opt(1)]';
MicA_MicVA = [TTopt(2),TTmax(2),TTmin(2),Mopt(2),TToptR0(2),TTR01(2),R0opt(2)]';
MicC_MicVC = [TTopt(3),TTmax(3),TTmin(3),Mopt(3),TToptR0(3),TTR01(3),R0opt(3)]';
CP = table(Cardinal_Parameters,MicB_MicVB,MicA_MicVA,MicC_MicVC)


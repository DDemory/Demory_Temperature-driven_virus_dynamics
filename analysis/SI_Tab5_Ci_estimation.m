%% A thermal trade-off between viral production and degradation drives phytoplankton-virus population dynamics
% === Code for: ===
% SI Tab 5 -- Estimation of the hyperparameter 95% confidence intervals
% David Demory - Jan 2021

% In this code we estimated the 95% confidence intervals by approximating
% the inverse of the Hessian matrix and calculating the variance matrix
% with the sandwich estimator method.
% Note: for parameters that are big or close to 0, the inverse matrix can
% be close to singular. These intervals need to be taken with caution.

%% ========== SETUP ==========
% (1) Which strain? 829, 451 or 834
clear all;close all;clc;

strain = 834;
disp(['Strain used = ',num2str(strain)])

% Parameters names as in the manuscript
paraname = {'A1','E1','A2','E2','K','\phi_K','T_\phi','\phi_r','s1','d1',...
    's2','d2','b1','T_\epsilon','\epsilon_r','\sigma_1','\delta_1',...
    '\delta_2','\omega'};

% option for the ode integration
options = odeset(); % Given your matlab version, you will need to modify
% the tolerance of the integration and the dp values. This should not
% significantly affect the confidence interval estimates for small
% enough epsilon. This code has been generated with Matlab2019b.


% load best parameters and data
dt = 0.1;
tend = 5;
time = [0:dt:tend]; % time for the ode integration;

% deltap (variation on the optimal parameters)
if strain == 829;
    load('pmin_829_Arrhenius_v2.mat')
    dp = 1E-11*pmin;

    % Treatment and variable names
    t_names = {'t_95C','t_125C','t_20C','t_25C','t_275C','t_30C'};
    V_names = {'V_95C','V_125C','V_20C','V_25C','V_275C','V_30C'};
    H_names = {'H_95C','H_125C','H_20C','H_25C','H_275C','H_30C'};
    tV_names = {'tV_95C','tV_125C','tV_20C','tV_25C','tV_275C','tV_30C'};
    tH_names = {'tH_95C','tH_125C','tH_20C','tH_25C','tH_275C','tH_30C'};
    y_names = {'y_95C','y_125C','y_20C','y_25C','y_275C','y_30C'};

elseif strain == 451;
    load('pmin_451_Arrhenius_v2.mat')
    dp = 5E-13*pmin;

    % Treatment and variable names
    t_names = {'t_125C','t_20C','t_25C','t_275C','t_30C'};
    V_names = {'V_125C','V_20C','V_25C','V_275C','V_30C'};
    H_names = {'H_125C','H_20C','H_25C','H_275C','H_30C'};
    tV_names = {'tV_125C','tV_20C','tV_25C','tV_275C','tV_30C'};
    tH_names = {'tH_125C','tH_20C','tH_25C','tH_275C','tH_30C'};
    y_names = {'y_125C','y_20C','y_25C','y_275C','y_30C'};

elseif strain == 834
    load('pmin_834_Arrhenius_v2.mat')
    dp = 1E-12*pmin;
    dp(1) = 5E-15*pmin(1);
    dp(2) = 3E-16*pmin(2);
    dp(5) = 4E-15*pmin(5);
    dp(17:18) = 1E-14*pmin(17:18);

    % Treatment and variable names
    t_names = {'t_125C','t_20C','t_25C','t_275C','t_30C'};
    V_names = {'V_125C','V_20C','V_25C','V_275C','V_30C'};
    H_names = {'H_125C','H_20C','H_25C','H_275C','H_30C'};
    tV_names = {'tV_125C','tV_20C','tV_25C','tV_275C','tV_30C'};
    tH_names = {'tH_125C','tH_20C','tH_25C','tH_275C','tH_30C'};
    y_names = {'y_125C','y_20C','y_25C','y_275C','y_30C'};

end

% hyperparameter id
i_param=[1:1:length(pmin)];
% load data
data = experiments(strain,dt,tend);


%% ========== Estimation of the sensitivity matrice ============

% storage variables
H_nom = []; V_nom = [];
Cvec_H = []; Cvec_V = [];
D = [];dd = [];
DHall = [];
DVall = [];
Difftot_V = [];
Difftot_H = [];

% temperature treatments related to the experiments in Demory et al. 2017
if strain == 829 % 6 treatments for 829 (9.5,12.5,20,25,27.5 and 30C)
    jend = 6;
else
    jend = 5; % 5 treatments for 451 and 834 (12.5,20,25,27.5 and 30C)
end

% ---------- loop on treatments ----------
for j = 1:jend

    % temperature of treatment j
    T = data.T(j);

    % ---------- CLEANING DATA ----------
    % looking for times were both H and V are measured
    tHi = data.(tH_names{j}); % time host data
    tVi = data.(tV_names{j}); % time virus data
    % Data
    H_data = data.(H_names{j});   % concentration host data
    V_data = data.(V_names{j});   % concentration virus data

    tHi(isnan(H_data)==1)=[];
    H_data(isnan(H_data)==1)=[];
    tVi(isnan(V_data)==1)=[];
    V_data(isnan(V_data)==1)=[];

    % find the common positive times in tH and tV
    indH=[];    indV=[];
    for iii=1:length(tHi),

        test= (abs(tVi-tHi(iii)) <1E-3);
        if  sum(test),

            if ~(tHi(iii)<0),
                indH=[indH;iii];

                indVi=find( test);
                indV=[indV;indVi];
            end
        end

    end

    tV=tVi(indV);
    tH=tHi(indH);

    yH = H_data(indH);
    yV = V_data(indV);


    % ---------- Calcul of the error's variance (Sig2) ----------
    % error of the fits using the optimal reference parameters (pmin)

    % initial conditions for the ode integration
    iniC = [yH(1),0,yV(1),0];

    % log-transform parameters
    p_ref = log(pmin);         % reference optinal parameters

    % integration
    [tmod,y_nominal] = ode15s(@SIV2_v2_log,time,iniC,options,p_ref,T);

    % results of the integration
    ysolH = y_nominal(:,1)+y_nominal(:,2); %Host (S+I)
    ysolV = y_nominal(:,3)+y_nominal(:,4); %Virus (Vi+Vni)

    % take only the value at time point similar for the host and virus data
    modH = interp1(tmod,ysolH,tH);
    modV = interp1(tmod,ysolV,tV);

    % ---------- Calcul for the sig2 ----------
    difH = (yH-modH).^2; % residual sum of squares (RSS) for the host
    difV = (yV-modV).^2; % RSS for the virus
    % Save all RSS (calculated for each temperature) in a matrix
    Difftot_H = [Difftot_H;difH];
    Difftot_V = [Difftot_V;difV];

    % ---------- Calcul the sensitivity matrice ----------
    % Here we compare the outputs from the optimal reference parameter with
    % parameters that vary around the optimal.

    % storage variables
    DH = []; DV =[];
    DH = []; DV =[];

    % nbr of parameter in the analysis
    nb_param=length(i_param);
    l_param=length(i_param);

    for i = i_param, %% loop on the parameters

        % change the parameter i by a factor deltap

        % sensitivity parameters (pmin+dp)
        p_sens    = log(pmin);
        p_sens(i) = log(pmin(i)+dp(i)); % variation around parameter i

        % storage variable
        sensH=[0];        sensV=[0];

        % Compute the sensitivity from each data points
        for ikp=2:length(tH)
            windowsensi=1;
            timesensi=[tH(ikp) tH(ikp)+windowsensi];

            % Integration using the reference pmin and p_sens
            [tmod,y_ref] = ode15s(@SIV2_v2_log,timesensi,y_nominal(ikp,:),options,p_ref,T);
            [tmod,y_sens] = ode15s(@SIV2_v2_log,timesensi,y_nominal(ikp,:),options,p_sens,T);

            % results of the integration
            ysolH_sens = y_sens(2,1)+y_sens(2,2); % new Host
            ysolV_sens = y_sens(2,3)+y_sens(2,4); % new Virus
            ysolH_ref = y_ref(2,1)+y_ref(2,2);    % reference Host
            ysolV_ref = y_ref(2,3)+y_ref(2,4);    % reference Virus

            % Calcul of the sensitivity (dy/dp)
            p = pmin(i);                  % parameter i
            dyH = (ysolH_sens-ysolH_ref); % difference host dH
            dyV = (ysolV_sens-ysolV_ref); % difference virus dV

            %dy/d(log(p))=((y(log(p+deltap)) - y(log(p)))/ deltap)*p
            sensH=[sensH;p*(dyH/dp(i))]; %dy/d(log(p)) = p*(dy/dp)
            sensV=[sensV;p*(dyV/dp(i))];
        end

        % store all the sensitivity in a vector for each treatments
        DH=[DH sensH];
        DV=[DV sensV];

    end

    % store the sensitivity into a matrix of sensitivity
    DHall=[DHall ; DH ];
    DVall=[DVall ; DV] ;
end

%% ========== Calcul of the Jacobian and the Hessian ==========

% Compute the Jacobian of f of log(p)
% Here I concatenate both sensitivity for the host and the virus to have
% one matrix and one Jacobian
J = [DHall;DVall]; % dim = n*p

% Compute the Hessian
H = (J'*J); % dim = p*p

%% ========== Compute the confidence intervals ==========

% Compute the diagonal of the RSS
% Same as the jacobian: I concatenate both Host and Virus RSS.
Difftot = [Difftot_H;Difftot_V];

% Compute the Variance Matrice Using the sandwich estimator
C = inv(J'*J)*J'*diag(Difftot)*J*inv(J'*J); % dim = p*p matrix
diagC = diag(C);

% Confidence intervals
% log Confidence interval
Cil = p_ref(i_param)' - 1.96*sqrt(diagC);
Ciu = p_ref(i_param)' + 1.96*sqrt(diagC);
ConfI_sandwitch_log = [Cil,Ciu];

% Confidence interval 95%
% We need to come back to the linear space: exp(logCI).
ConfI_sandwitch_p = exp(ConfI_sandwitch_log);


%% Show 95% Ci results

parameter = paraname(i_param)';

log_low = ConfI_sandwitch_log(:,1);
log_optimal = p_ref(i_param)';
log_up = ConfI_sandwitch_log(:,2);
diagC  = diag(C);
Std = exp(sqrt(diag(C)));
diagH  = diag(H);
id     = i_param';
CIlog = table(id,parameter,diagC,log_low, log_optimal, log_up);

low = ConfI_sandwitch_p(:,1);
optimal = pmin(i_param)';
up = ConfI_sandwitch_p(:,2);

CI95 = table(id,parameter,low, optimal, up)

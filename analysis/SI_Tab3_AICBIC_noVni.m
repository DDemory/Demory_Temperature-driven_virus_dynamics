%% A thermal trade-off between viral production and degradation drives phytoplankton-virus population dynamics
% === Code for: ===
% Calculation of AIC and BIC for the temperature-driven model and the basal
% model for the pair MicB/MicV-B.
% David Demory -- Jan 2021

%% Load parameters and data
strain = 829; % Choose 829 (MicB/MicV-B)

if strain ~= 829
    disp('Wrong strain number, choose strain = 829, 451 or 834')
    return
end

d0 = 0;
df = 5;
dt = 0.01;

% Load experimental data
y_data = experiments(strain,dt,d0);

% Number of treatments
nexp = 6;

% Treatment and variable names
t_names = {'t_95C','t_125C','t_20C','t_25C','t_275C','t_30C'};
V_names = {'V_95C','V_125C','V_20C','V_25C','V_275C','V_30C'};
H_names = {'H_95C','H_125C','H_20C','H_25C','H_275C','H_30C'};
tV_names = {'tV_95C','tV_125C','tV_20C','tV_25C','tV_275C','tV_30C'};
tH_names = {'tH_95C','tH_125C','tH_20C','tH_25C','tH_275C','tH_30C'};
sdV_names = {'sdV_95C','sdV_125C','sdV_20C','sdV_25C','sdV_275C','sdV_30C'};
sdH_names = {'sdH_95C','sdH_125C','sdH_20C','sdH_25C','sdH_275C','sdH_30C'};
y_names = {'y_95C','y_125C','y_20C','y_25C','y_275C','y_30C'};

% Experimental temperatures
TT = [9.5,12.5,20,25,27.5,30];

% Model parameters
load('pmin_829_Arrhenius_v2.mat')

% pmin linear decay
%pmin =[1.1334e+10, 6743, 1.2515e+25, 17349, 1.0013e+09, 1.0808e-06, 39.676, 0.14913, 1.0126e+22,...
%    13119, 1.6263e+22, 13264, 17.021, 297.62, 6.5593, 1.2344e+19, 0.001863, 0, ...
%    2.0012e-08];

% Load noVni model parameters
load('pmin_829_noVni_v2.mat')

% time of the simulation
t = d0:dt:df;

sumSSE=0;

hostfigid = [1,3,5,7,9,11];
virusfigid = [2,4,6,8,10,12];

% color
lightgrey = [0.85 0.85 0.85];
colH      = [26, 188, 156 ]/255;
colV      = [165, 105, 189]/255;

for i = 1:nexp;
    
    % Temperature of the treament i
    T = TT(i);
    
    % Data
    H_data = y_data.(H_names{i});   % concentration host data
    tH_data = y_data.(tH_names{i}); % time host data
    V_data = y_data.(V_names{i});   % concentration virus data
    tV_data = y_data.(tV_names{i}); % time virus data
    
    tH_data(isnan(H_data)==1)=[];
    H_data(isnan(H_data)==1)=[];
    tV_data(isnan(V_data)==1)=[];
    V_data(isnan(V_data)==1)=[];
    
    % ODE Integration
    % initial condition
    Hinistat_std=std(H_data(1:3));   % std initial condition host
    Vinistat_std=std(V_data(1:3));   % std initial condition virus
    Hinistat_mean=mean(H_data(1:3)); % mean initial condition host
    Vinistat_mean=mean(V_data(1:3)); % mean initial condition virus
    
    % ODE options
    options = odeset('reltol',1E-3,'abstol',1E-3);
    
    % Variable for the fill areas
    yHplus=[]; yHmin = [];yVmoins=[]; yVplus = [];
    
    % Temperature variation \pm 1
    dx = 0.001; % paper version
    dx = 0.1;
    
    for deltaT=-1:dx:1
        
        % host initial condition
        Hini=max(0,Hinistat_mean+randn(1)*Hinistat_std);
        
        % virus initial condition
        Vini=max(0,Vinistat_mean+0.3*randn(1)*Vinistat_std);
        
        % if no variation of temperature the initial condition is equal to
        % the experimental concentration at t0
        if deltaT==0
            Hini=H_data(1);
            Vini=V_data(1);
        end
        
        % initial conditions for ode45
        Ci = [Hini,0,Vini,0];
        
        % Calculate the parameters for T+deltaT
        [mu,K,phi,lambda,epsilon,beta,sigma,delta,omega,m]=Tdriven_fct(T+deltaT,pmin);
        
        % Integration using ode45 at T+deltaT
        [tmod,ymod] = ode45(@SIV2_v2,t,Ci,options,pmin,T+deltaT);
        
        % Use the H0 for the comparison with strain 829
        if strain == 829
            % Integration
            [tfit0,yfit0] = ode45(@SIV_v2_noVni,t, [H_data(1),0,V_data(1)], [], pmin_noVni, T);
            % Total Host = S+I
            ysol0H = yfit0(:,1)+yfit0(:,2);
            % Total Virus = Vi+Vni
            ysol0V = yfit0(:,3);
        end
        
        % Calcul on the envelops
        if deltaT==0;
            simnominaleH=ymod(:,1)+ymod(:,2);
            simnominaleV=ymod(:,3)+ymod(:,4);
        end
    end
    
    %% AIC and BIC estimations
    
    %% AIC and BIC estimations
    
    % Interpolation ysim
    ysim_intH=interp1(tmod,simnominaleH,tH_data);
    ysim_intV=interp1(tmod,simnominaleV,tV_data);
    ysim_int0H=interp1(tmod,yfit0(:,1)+yfit0(:,2),tH_data);
    ysim_int0V=interp1(tmod,yfit0(:,3),tV_data);
    
    % SSE
    ecartR_H = sum((ysim_intH-H_data).^2);
    ecartR_V = sum((ysim_intV-V_data).^2);
    ecartR0_H = sum((ysim_int0H-H_data).^2);
    ecartR0_V = sum((ysim_int0V-V_data).^2);
    
    % number of parameters
    %k0 = 10; % Basal model
    k0 = 16; % noVni model
    kT = 19; % Temperature-driven model
    
    % data size
    nH = length(H_data);
    nV = length(V_data);
    
    % AIC
    AIC_T = 2*kT + nH*log(ecartR_H/nH) + nV*log(ecartR_V/nV);
    AIC_0 = 2*k0 + nH*log(ecartR0_H/nH) + nV*log(ecartR0_V/nV);
    
    % BIC
    BIC_T = kT*log(nH+nV) + nH*log(ecartR_H/nH) + nV*log(ecartR_V/nV);
    BIC_0 = k0*log(nH+nV) + nH*log(ecartR0_H/nH) + nV*log(ecartR0_V/nV);
    
    % Print AIC and BIC results
    disp('################################################################################')
    disp(sprintf(' Temperature:  %1.2f',T))
    disp(sprintf(' Temperature-driven model:  AIC_T = %0.5g,   BIC_T =%0.5g',[AIC_T BIC_T]))
    disp(sprintf(' no_Vni model:  AIC_0 = %0.5g,   BIC_0 =%0.5g',[AIC_0 BIC_0]))
end



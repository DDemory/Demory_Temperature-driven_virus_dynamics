function [data,Ci,tspan,T] = experiments(Mstrain,dt,tend)
% Load experimental data
% David Demory
% Mstrain: Micromonas-virus pairs: 829 (Mic-B/MicV-B), 834 (Mic-C/MicV-C) or 451 (Mic-A/MicV-A)
% dt = time step for integration
% tend = final time for integration

%% ===== Data folder =====
data_folder = 'data/';

if Mstrain == 829
    fileV = 'data_B829_VIRUS.xlsx';
    fileH = 'data_B829_HOST.xlsx';
    nexp = 6;
    t_names = {'t_95C','t_125C','t_20C','t_25C','t_275C','t_30C'};
    V_names = {'V_95C','V_125C','V_20C','V_25C','V_275C','V_30C'};
    H_names = {'H_95C','H_125C','H_20C','H_25C','H_275C','H_30C'};
    tV_names = {'tV_95C','tV_125C','tV_20C','tV_25C','tV_275C','tV_30C'};
    tH_names = {'tH_95C','tH_125C','tH_20C','tH_25C','tH_275C','tH_30C'};
    T = [9.5,12.5,20,25,27.5,30];

elseif Mstrain == 451
    fileV = 'data_A451_VIRUS.xlsx';
    fileH = 'data_A451_HOST.xlsx';
    nexp = 5;
    t_names = {'t_125C','t_20C','t_25C','t_275C','t_30C'};
    V_names = {'V_125C','V_20C','V_25C','V_275C','V_30C'};
    H_names = {'H_125C','H_20C','H_25C','H_275C','H_30C'};
    tV_names = {'tV_125C','tV_20C','tV_25C','tV_275C','tV_30C'};
    tH_names = {'tH_125C','tH_20C','tH_25C','tH_275C','tH_30C'};
    T = [12.5,20,25,27.5,30];

elseif Mstrain == 834
    fileV = 'data_C834_VIRUS.xlsx';
    fileH = 'data_C834_HOST.xlsx';
    nexp = 5;
    t_names = {'t_125C','t_20C','t_25C','t_275C','t_30C'};
    V_names = {'V_125C','V_20C','V_25C','V_275C','V_30C'};
    H_names = {'H_125C','H_20C','H_25C','H_275C','H_30C'};
    tV_names = {'tV_125C','tV_20C','tV_25C','tV_275C','tV_30C'};
    tH_names = {'tH_125C','tH_20C','tH_25C','tH_275C','tH_30C'};
    T = [12.5,20,25,27.5,30];

end


for i = 1:nexp;

    dataV = xlsread([data_folder fileV],i);
    dataH = xlsread([data_folder fileH],i); 

    tV = dataV(:,1);
    V = dataV(:,2);

    tH = dataH(:,1);
    H = dataH(:,2);

    data.(V_names{i}) = V;
    data.(H_names{i}) = H;
    data.(tV_names{i}) = tV;
    data.(tH_names{i}) = tH;
    data.T = T;

    Ci.(V_names{i}) = V(1);
    Ci.(H_names{i}) = H(1);

    t = 0:dt:tend;
    tspan.(t_names{i}) =t;
end
end

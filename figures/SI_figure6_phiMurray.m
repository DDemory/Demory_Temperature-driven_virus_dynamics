%% A thermal trade-off between viral production and degradation drives phytoplankton-virus population dynamics
% === Code for: ===
% SI Figure 4 -- Biophysical maximum adsorption rate based on Murray and Jackson 1992
% David Demory - Jan 2021

%% Calculus of Murray and Jackson maximum biophysical adsorption.

% Température
T = 1:1:30;    % Celcius
TK = 273.15+T; % Kelvin

% Viscosity
dynVisco = [0.001843,0.001783,0.001726,0.001671,0.001620,0.001571,0.001524,...
    0.001480,0.001438,0.001397,0.001359,0.001322,0.001286,0.001252,0.001220,0.001189,...
    0.001159,0.001131,0.001103,0.001077,0.001051,0.001027,0.001004,0.000981,...
    0.000959,0.000938,0.000918,0.000898,0.000879,0.000861]; %\mu viscosity [ kg/m-s]
eta = dynVisco.*1e-6;

% Botzman constant
Kb = 1.380649*10^-23; % Boltzmann constant [ kg⋅m^2⋅s^−2⋅K^−1 ]
Kb = Kb * 10^12; % Boltzmann constant [ kg⋅micron^2⋅s^−2⋅K^−1 ]

% Radius
d_Micro = 1.5; % μm
r_Micro = d_Micro/2;
d_Virus = 130 * 0.001; % μm
r_Virus = d_Virus/2;

% Einstein formulations for host and viruses
DMicro = (Kb.*TK)./(6*pi*eta*r_Micro);
DVirus = (Kb.*TK)./(6*pi*eta*r_Virus);

% Adsorption
a1 = (4*pi *(DMicro + DVirus).*(r_Micro + r_Virus)); %[um^3 s^-1]
phi_max = a1*(10^-12)*86400; %[mL day^-1]

%% Plot
hfig = figure('color','White');
plot(T,phi_max,'LineWidth',3);
xlabel('Temperature')
ylabel('Maximum adsorption (mL/d)')
xlim([0 31])
set(gca,'FontSize',14)
set(gcf,'Color','white')


% Fontsize
set(findall(gcf,'-property','FontSize'),'FontSize',18)
% Save
print(hfig,'-depsc','-r600','outputs/SI_Figure6_phiMurray.eps')

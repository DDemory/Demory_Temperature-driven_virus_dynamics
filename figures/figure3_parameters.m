%% A thermal trade-off between viral production and degradation drives phytoplankton-virus population dynamics
% Figure 3 -- Temperature-driven functions for MicB/MicV-B (strain 829)
% David Demory

%% Temperature-driven functions
% Temperature range
T = [0:0.01:32.5];
% Hyperparameters of strain 829
load('pmin_829.mat')
% Calculate model parameters at T
[mu,K,phi,lambda,epsilon,beta,sigma,delta,omega,m] = Tdriven_fct(T,pmin);

%% DATA
% Host growth from Demory et al. ISME 2018.
% Experimental temperature
Tmu = [4,4,7.5,7.5,7.5,7.5,7.5,7.5,12.5,12.5,12.5,12.5,12.5,12.5,15,15,20,20,20,20,20,25,25,25,25,25,25,25,27.5,27.5,27.5,27.5,27.5,27.5,30,30,30,30,32.5,32.5,32.5,32.5,32.5];
% Host net growth rate
muexp = [0,0.0325,0,0,0,0.135,0.2705,0,0.2415,0.2442,0.3189,0.3922,0.2876,0.4129,0.48,0.6014,0.832,0.9168,0.7127,0.8344,0.8699,1.0091,1.1284,1.1773,0.9358,0.9799,1.0392,...
    1.1228,0.9808,0.9833,1.2815,1.0546,1.0598,1.1374,0.9743,0.9932,1.0049,1.0142,0,0,0,0,0];

% Virus parameters from Demory et al. ISME 2017.
% Latent period temperatures
Texp = [9.5,12.5,20,25,27.5,30];
% Latent period max
Lexp = ([11,11,7,3,3,7]);
% Latent period min
Lexp2 = ([7,7,3,1,1,3]);
% stats
L = (Lexp+Lexp2)/2;
sL = std([Lexp;Lexp2]);
% Burst size
Bexp = [NaN,84,142,139,177,49];
% host growth
muexp2 = [0.5937,0.6268,0.8909,1.0895,1.1152,0.8721];
% decay temperatures
Texpd = [4,12.5,20,25,27.5,30];
% decay rate
dexps = [0.0845,0.1155,0.1784,0.1652,0.2385,0.2344];
% lost of infectivity
dexpi = [0.1511,0.2045,0.2903,0.2839,0.3937,0.3928];

%% Plot
hfig = figure('position',[0 0 1100 600]);
hfig.Color = 'white';

% Host growth rate \mu
subplot(2,3,1)
hold on
plot(Tmu,muexp,'.','Markersize',40,'Color',[0.8 0.8 0.8])
plot(Texp,muexp2,'k.','Markersize',40)
plot(T,mu-m,'k-','linewidth',3,'linestyle','-')
ylabel('day^{-1}')
set(gca,'Xlim',[0 35],'xtick',[0,10,20,30])
set(gca,'FontSize',18,'FontName', 'Times New Roman')
set(gca,'TickDir','out')
ax = gca;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 18;
ax.TickDir = 'out';
title('Net growth rate ($\mu - \psi$)','Interpreter','latex')

% Latent period \tau
subplot(2,3,2)
hold on
plot(Texp,L/24,'k.','Markersize',40)
text(0.02,2.5,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','FontSize',28)
plot(T,(1./lambda),'k-','linewidth',3,'linestyle','-')
ylabel('hours')
set(gca,'Xlim',[0 35],'xtick',[0,10,20,30])
set(gca,'Ylim',[0 0.7],'ytick',[0,0.21,0.42,0.63])
set(gca,'FontSize',18,'FontName', 'Times New Roman')
set(gca,'TickDir','out')
ax = gca;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 18;
ax.TickDir = 'out';
title('Latente period ($\tau$)','Interpreter','latex')

% Burst size \beta
subplot(2,3,3)
hold on
plot(Texp,Bexp,'k.','Markersize',40)
plot(T,beta,'k-','linewidth',3,'linestyle','-')
ylabel('Virus/Lysed cell')
set(gca,'Xlim',[0 35],'xtick',[0,10,20,30])
set(gca,'Ylim',[0 200],'ytick',[0,65,130,195])
set(gca,'FontSize',18,'FontName', 'Times New Roman')
set(gca,'TickDir','out')
ax = gca;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 18;
ax.TickDir = 'out';
title('Viral burst size ($\beta$)','Interpreter','latex')

% Fraction of non-infectious virus \epsilon
subplot(2,3,4)
plot(T,epsilon,'k-','linewidth',3,'linestyle','-')
ylabel('Fraction')
set(gca,'Xlim',[0 35],'xtick',[0,10,20,30])
set(gca,'Ylim',[0 1],'ytick',[0,0.3,0.6,0.9])
set(gca,'FontSize',18,'FontName', 'Times New Roman')
set(gca,'TickDir','out')
ax = gca;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 18;
ax.TickDir = 'out';
title('Fraction of non-infectious virus ($\epsilon$)','Interpreter','latex')

% Decay rates \delta and \sigma
subplot(2,3,5)
hold on
plot(Texpd,dexps,'k.','Markersize',40)
plot(Texpd,dexpi,'k^','Markersize',10)
plot(T,delta,'k-','linewidth',3,'linestyle','-')
plot(T,sigma,'k--','linewidth',3,'linestyle','--')
ylabel('day^{-1}')
set(gca,'Xlim',[0 35],'xtick',[0,10,20,30])
set(gca,'Ylim',[0 2],'ytick',[0,0.7,1.4,2])
set(gca,'FontSize',18,'FontName', 'Times New Roman')
set(gca,'TickDir','out')
ax = gca;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 18;
ax.TickDir = 'out';
title('Decay rate($\delta$ and $\sigma$)','Interpreter','latex')
xlabel('Temperature (^{\circ}C)')

% Adsorption rate \phi
subplot(2,3,6)
plot(T,phi,'k-','linewidth',3,'linestyle','-')
ylabel('mL cell^{-1} day^{-1}')
set(gca,'Xlim',[0 35],'xtick',[0,10,20,30])
set(gca,'Ylim',[0 0.45E-6],'ytick',[0,0.15,0.3,0.45]*1E-6)
set(gca,'FontSize',18,'FontName', 'Times New Roman')
set(gca,'TickDir','out')
ax = gca;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.FontSize = 18;
ax.TickDir = 'out';
title('Adsorption rate ($\phi$)','Interpreter','latex')

% Save figure
print(hfig,'-depsc','-r600','outputs/Figure3_parameters.eps')


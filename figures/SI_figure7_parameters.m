%% A thermal trade-off between viral production and degradation drives phytoplankton-virus population dynamics
% === Code for: ===
% SI Figure 5 -- Life history trait temperature functions for the 3 pairs
% David Demory - Jan 2021

%% setup

% Temperature
T = 0:0.1:35;
    
% load parameters
load('pmin_829_Arrhenius_v2.mat')
[mu3,K3,phi3,lambda3,epsilon3,beta3,sigma3,delta3,omega3,psi3]=Tdriven_fct(T,pmin);
clear pmin
load('pmin_451_Arrhenius_v2.mat')
[mu1,K1,phi1,lambda1,epsilon1,beta1,sigma1,delta1,omega1,psi1]=Tdriven_fct(T,pmin);
clear pmin
load('pmin_834_Arrhenius_v2.mat')
[mu2,K2,phi2,lambda2,epsilon2,beta2,sigma2,delta2,omega2,psi2]=Tdriven_fct(T,pmin);
clear pmin

%% Plot
hfig = figure('color','White','position',[0 0 1000 500]);

% mu and psi
subplot(2,3,1)
hold on
B = plot(T,mu3-psi3,'LineWidth',3);
cB = get(B,'Color');
A = plot(T,mu1-psi1,'LineWidth',3);
cA = get(A,'Color');
C = plot(T,mu2-psi2,'LineWidth',3);
cC = get(C,'Color');
ylim([0 1.1])
xlim([0 40])
legend('Mic-B/MicV-B','Mic-A/MicV-A','Mic-C/MicV-C','location','northoutside')
legend boxoff
title('Net Growth rate (\mu - \psi)')
ylabel('d^{-1}')

% lambda
subplot(2,3,2)
hold on
plot(T,lambda3,'LineWidth',3)
plot(T,lambda1,'LineWidth',3)
plot(T,lambda2,'LineWidth',3)
title('Lysis rate (\lambda)')
ylabel('d^{-1}')

% beta
subplot(2,3,3)
hold on
plot(T,beta3,'LineWidth',3)
plot(T,beta1,'LineWidth',3)
plot(T,beta2,'LineWidth',3)
title('Burst size (\beta)')
ylabel('Virus/cell')

% epsilon
subplot(2,3,4)
hold on
plot(T,epsilon3,'LineWidth',3)
plot(T,epsilon1,'LineWidth',3)
plot(T,epsilon2,'LineWidth',3)
title('% of non-infectious virus (\epsilon)')


subplot(2,3,5)
hold on
%sigma
plot(T,sigma3,'LineWidth',3)
plot(T,sigma1,'LineWidth',3)
plot(T,sigma2,'LineWidth',3)
a = plot(0,0,'k-','LineWidth',3);
%delta
plot(T,delta3,'--','LineWidth',3,'color',cB)
plot(T,delta1,'--','LineWidth',3,'color',cA)
plot(T,delta2,'--','LineWidth',3,'color',cC)
b = plot(0,0,'k--','LineWidth',3);
set(gca,'Yscale','log')
legend([a,b],{'\sigma','\delta'},'location','best')
legend boxoff
title('Decay rates (\delta and \sigma)')
ylabel('log d^{-1}')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
xlabel('Temperature (^\circC)')

subplot(2,3,6)
hold on
plot(T,phi3,'LineWidth',3)
plot(T,phi1,'LineWidth',3)
plot(T,phi2,'LineWidth',3)
set(gca,'Yscale','log')
title('Adsorption rate (\phi)')
ylabel('mL cell^{-1} d^{-1}')

% Save
print(hfig,'-depsc','-r600','outputs/SI_Figure7_parameters.eps')

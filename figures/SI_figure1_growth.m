%% A thermal trade-off between viral production and degradation drives phytoplankton-virus population dynamics
% === Code for: ===
% SI Figure 9 -- growth and mortality function
% David Demory - Jan 2021

%% Setup
T = 0:0.1:35;

load('pmin_829_Arrhenius_v2.mat')
[mu1,K1,phi1,lambda1,epsilon1,beta1,sigma1,delta1,omega1,m1]=Tdriven_fct(T,pmin);
clear pmin

load('pmin_451_Arrhenius_v2.mat')
[mu2,K2,phi2,lambda2,epsilon2,beta2,sigma2,delta2,omega2,m2]=Tdriven_fct(T,pmin);

load('pmin_834_Arrhenius_v2.mat')
[mu3,K3,phi3,lambda3,epsilon3,beta3,sigma3,delta3,omega3,m3]=Tdriven_fct(T,pmin);


%% Plot

hfig = figure('color','White','position',[0 0 1000 250]);

% Net growth
subplot(1,3,1)
hold on
A = plot(T,mu1-m1,'LineWidth',3);
cA = get(A,'Color');
B = plot(T,mu2-m2,'LineWidth',3);
cB = get(B,'Color');
C = plot(T,mu3-m3,'LineWidth',3);
cC = get(C,'Color');
ylim([0 1.1])
xlim([0 40])
legend('Mic-B','Mic-A','Mic-C','location','NorthOutside')
legend boxoff
title('Net Growth')
ylabel('d^{-1}')

% Gross growth
subplot(1,3,2)
hold on
A = plot(T,mu1,'LineWidth',3);
cA = get(A,'Color');
B = plot(T,mu2,'LineWidth',3);
cB = get(B,'Color');
C = plot(T,mu3,'LineWidth',3);
cC = get(C,'Color');
%ylim([0 1])
xlim([0 40])
title('Gross Growth')
xlabel('Temperature (^\circ C)')

% psi
subplot(1,3,3)
hold on
A = plot(T,m1,'LineWidth',3);
cA = get(A,'Color');
B = plot(T,m2,'LineWidth',3);
cB = get(B,'Color');
C = plot(T,m3,'LineWidth',3);
cC = get(C,'Color');
xlim([0 40])
title('Non-lysis loss')

set(findall(gcf,'-property','FontSize'),'FontSize',12)

% Save
print(hfig,'-depsc','-r600','outputs/SI_Figure1_growth.eps')

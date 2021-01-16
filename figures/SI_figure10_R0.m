%% A thermal trade-off between viral production and degradation drives phytoplankton-virus population dynamics
% === Code for: ===
% SI Figure 8 -- R0 as function of the temperature
% David Demory - Jan 2021

%% setup
T = [0:0.0001:33.5];
hfig = figure('position',[0 0 1000 500])
names = {'Mic-B/MicV-B','Mic-A/MicV-A','Mic-C/MicV-C'}
EpiNiche = []; RfgNiche = []; DiffNiche = []; R0max = []; TTopt = []; Mopt =[];

%% Temperature-driven functions
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


%% Calcul R0;
R0 = ((1-epsilon).*beta.*lambda.*phi.*K)./((lambda+psi).*(phi*K+sigma));

%% Plot
hold on
plot(T,R0,'LineWidth',3)
drawnow
end

%% save figure
legend('Mic-B/MicV-B','Mic-A/MicV-A','Mic-C/MicV-C','location','best')
legend boxoff
xlabel('Temperature $^\circ$C','Interpreter','latex')
ylabel('${\cal{R}}_0$','Interpreter','latex')
set(gca,'FontSize',16)
set(gcf,'color','White','position',[0 0 600 400])
% Save
print(hfig,'-depsc','-r600','outputs/SI_Figure10_R0.eps')

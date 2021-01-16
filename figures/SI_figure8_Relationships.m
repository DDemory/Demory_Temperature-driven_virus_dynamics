%% A thermal trade-off between viral production and degradation drives phytoplankton-virus population dynamics
% === Code for: ===
% SI Figure 6 -- Relationships between temperature parameters describing
% the thermal niche.
% David Demory - Jan 2021

%% setup

T = [0:0.0001:33.5];
names = {'Mic-B/MicV-B','Mic-A/MicV-A','Mic-C/MicV-C'};
EpiNiche = []; RfgNiche = []; DiffNiche = []; R0max = []; TTopt = []; Mopt =[];

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
    
    %% Calcul R0;
    
    R0 = ((1-epsilon).*beta.*lambda.*phi.*K)./((lambda+psi).*(phi*K+sigma));
    ToptR0 = T(find(R0 == max(R0)));
    tempT = T(find(round(R0,2)==1));TR01 = tempT(1);

    %% Calcul of temperature parameters
    %Epi zone
    EpiNiche = [EpiNiche,TR01];
    %Refuge zone
    RfgNiche = [RfgNiche,Tmax-TR01];
    % Diff
    DiffNiche = [DiffNiche,Topt-ToptR0];
    %R0max
    R0max = [R0max,max(R0)];
    %TTopt
    TTopt = [TTopt,Topt];
    % Muopt
    Mopt = [Mopt,muopt];
end

%% relationships
relation = figure('position',[0 0 500 1000],'Color','White');
subplot(3,1,1)
mdl = fitlm(DiffNiche,R0max)
plot(mdl)
legend off
title('')
xlabel('Optimum difference ($T_{opt}$ - $T_{R0,opt}$)','Interpreter','latex')
ylabel('$R0_{max}$','Interpreter','Latex')
R2adj = 0.001*round(mdl.Rsquared.Adjusted*1000);
pval = 0.001*round(mdl.anova.pValue(1)*1000);
text(1,10,['R$^2_{adj}$ = ',num2str(R2adj),' -- $p$-value = ',num2str(pval)],'Interpreter','Latex','Fontsize',12)
set(gca,'FontSize',16)
legend('Data','Fit','Condidence bounds','Fontsize',12,'location','Northoutside','orientation','horizontal')
legend boxoff
xlim([0 12])

subplot(3,1,2)
mdl = fitlm(DiffNiche,EpiNiche)
plot(mdl)
legend off
title('')
xlabel('Optimum difference ($T_{opt}$ - $T_{R0,opt}$)','Interpreter','Latex')
ylabel('Endemic niche','Interpreter','Latex')
R2adj = 0.001*round(mdl.Rsquared.Adjusted*1000);
pval = 0.001*round(mdl.anova.pValue(1)*1000);
text(1,24.5,['R$^2_{adj}$ = ',num2str(R2adj),' -- $p$-value = ',num2str(pval)],'Interpreter','Latex','Fontsize',12)
set(gca,'FontSize',16)
legend off
xlim([0 12])

subplot(3,1,3)
mdl = fitlm(DiffNiche,RfgNiche)
plot(mdl)
legend off
title('')
xlabel('Optimum difference ($T_{opt}$ - $T_{R0,opt}$)','Interpreter','Latex')
ylabel('Refuge niche','Interpreter','Latex')
R2adj = 0.001*round(mdl.Rsquared.Adjusted*1000);
pval = 0.001*round(mdl.anova.pValue(1)*1000);
text(5,5.5,['R$^2_{adj}$ = ',num2str(R2adj),' -- $p$-value = ',num2str(pval)],'Interpreter','Latex','Fontsize',12)
set(gca,'FontSize',16)
legend off
xlim([0 12])


%% Figure

hfig = figure('position',[0 0 500 1000],'Color','White');

subplot(3,1,1)
plot(DiffNiche,R0max,'ko','MarkerFaceColor','k','MarkerSize',10)
xlabel('Optimum difference ($T_{opt}$ - $T_{R0,opt}$)','Interpreter','latex')
ylabel('$R0_{max}$','Interpreter','Latex')
set(gca,'FontSize',16)
text(-2,170,'a','FontSize',18)
box off
xlim([0 12])

subplot(3,1,2)
plot(DiffNiche,EpiNiche,'ko','MarkerFaceColor','k','MarkerSize',10)
xlabel('Optimum difference ($T_{opt}$ - $T_{R0,opt}$)','Interpreter','Latex')
ylabel('Endemic niche','Interpreter','Latex')
set(gca,'FontSize',16)
text(-2,29,'b','FontSize',18)
box off
xlim([0 12])

subplot(3,1,3)
plot(DiffNiche,RfgNiche,'ko','MarkerFaceColor','k','MarkerSize',10)
xlabel('Optimum difference ($T_{opt}$ - $T_{R0,opt}$)','Interpreter','Latex')
ylabel('Refuge niche','Interpreter','Latex')
set(gca,'FontSize',16)
text(0.65,11,'c','FontSize',18)
box off

% Save
print(hfig,'-depsc','-r600','outputs/SI_Figure8_relationships.eps')





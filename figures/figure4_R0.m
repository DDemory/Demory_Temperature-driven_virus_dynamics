%% A thermal trade-off between viral production and degradation drives phytoplankton-virus population dynamics
% Figure 4 -- Ecological states as function of temperature (only fig. b and
% c.)
% David Demory

%% FONCTION DE LA TEMPERATURE -------------------------------------------------
% Temperature range
T = [0:0.0001:33.5];

% Title names
names = {'Mic-B/MicV-B','Mic-A/MicV-A','Mic-C/MicV-C'};

% Thermal niche initiation
EpiNiche = []; RfgNiche = []; DiffNiche = []; R0max = []; TTopt = []; Mopt =[];

% figure initiaiton 
hfig = figure('color','white','position',[0 0 1000 500]);

% loop on the strains
for s = 1:3
    % Which strain?
    % MicB/MicV-B
    if s == 1
        load('pmin_829.mat')
        [mu,K,phi,lambda,epsilon,beta,sigma,delta,omega,psi] = Tdriven_fct(T,pmin);
        [Topt,Tmax,muopt,Tmin] = Calculate_CT(pmin);
    % MicA/MicV-A
    elseif s == 2
        load('pmin_451.mat')
        [mu,K,phi,lambda,epsilon,beta,sigma,delta,omega,psi] = Tdriven_fct(T,pmin);
        [Topt,Tmax,muopt,Tmin] = Calculate_CT(pmin);
    % MicC/MicV-C
    elseif s == 3
        load('pmin_834.mat')
        [mu,K,phi,lambda,epsilon,beta,sigma,delta,omega,psi] = Tdriven_fct(T,pmin);
        [Topt,Tmax,muopt,Tmin] = Calculate_CT(pmin);
    end


%% Calcul R0 and thermal niches
% R0 at T
R0 = ((1-epsilon).*beta.*lambda.*phi.*K)./((lambda+psi).*(phi*K+sigma));

% Optimal R0 (T for R0 = max(R0))
ToptR0 = T(find(R0 == max(R0)));
% Maximal temperature of the infection niche (T for R0 = 1)
tempT = T(find(round(R0,2)==1));TR01 = tempT(1);

%Endemic niche
EpiNiche = [EpiNiche,TR01];
%Refuge niche
RfgNiche = [RfgNiche,Tmax-TR01];
% Difference between host and virus niches
DiffNiche = [DiffNiche,Topt-ToptR0];
%R0 max
R0max = [R0max,max(R0)];
%Optimal niche temperatures
TTopt = [TTopt,Topt];
%Optimal growth
Mopt = [Mopt,muopt];

%% Thermal niches

subplot(2,3,s)
yyaxis left
hold on

% -- Niches
pe = patch([0 TR01 TR01 0], [1.1 1.1 0 0],'g'); % endemic niche
pr = patch([TR01 Tmax Tmax TR01], [1.1 1.1 0 0], 'g'); % refuge niche
pc = patch([Tmax max(xlim) max(xlim) Tmax], [1.1 1.1 0 0], 'r'); % habitat loss niche
pe.FaceColor = [135,206,250]/255;pe.FaceAlpha = 0.2;pe.EdgeColor = 'none';
pr.FaceColor = [60,179,113]/255;pr.FaceAlpha = 0.2;pr.EdgeColor = 'none';
pc.FaceColor = [255 102 102]/255;pc.FaceAlpha = 0.2;pc.EdgeColor = 'none';
yline(1.2,'k--','LineWidth',2)
set(gca,'box','off')

% -- Net growth rate
plot(T,mu-psi,'-','LineWidth',3)
% Optimal growth 
plot([Topt, Topt],[0,max(mu-psi)],'k--')
plot([Topt Topt],[max(mu-psi) max(mu-psi)],'k.','MarkerSize',20)
% Maximal growth
plot([Tmax Tmax],[0 0],'k.','MarkerSize',20)
set(gca,'Ycolor','black')

% Left yaxis label
if s == 1
ylabel('Net growth (d^{-1})')
end
% Left axis lims
ylim([0 1.1])
yticks([0,0.5,1])

% -- R0 
yyaxis right
% R0 
plot(T,R0,'-','LineWidth',3)
% Optimal R0
plot([ToptR0(1) ToptR0(1)],[0 ,max(R0)],'k--')
plot([ToptR0 ToptR0],[max(R0) max(R0)],'k.','MarkerSize',20)
% T for R0 = 1
plot([TR01 TR01],[1.1 1.1],'k.','MarkerSize',20)
plot([TR01(1) TR01(1)],[0 1],'k--')
set(gca,'Ycolor','black')
% tickdir
set(gca,'TickDir','out')

% Right Yaxis label
if s == 2
xlabel('Temperature (^{\circ}C)')
end
% X axis label
if s == 3
ylabel('\it{R0}')
end

% lims and ticks
ylim([0 200])
yticks([0,100,200])

% titles
title(names(s))

end

%% Thermal niche differences
% Epidemic niche
subplot(2,3,4)
bar(EpiNiche,'FaceColor','k','EdgeColor','none','FaceAlpha',0.5)
box off
ylim([0 30])
ylabel('Temperature (^{\circ}C)')
xticklabels({'Mic-B/MicV-B','Mic-A/MicV-A','Mic-C/MicV-C'})
xtickangle(45)
title('Epidemic Niche')
text(0.65,30,num2str(round(EpiNiche(1),1)))
text(1.75,25,num2str(round(EpiNiche(2),1)))
text(2.65,27,num2str(round(EpiNiche(3),1)))

% Refuge Niche
subplot(2,3,5)
bar(RfgNiche,'FaceColor','k','EdgeColor','none','FaceAlpha',0.5)
box off
ylim([0 30])
xticklabels({'Mic-B/MicV-B','Mic-A/MicV-A','Mic-C/MicV-C'})
xtickangle(45)
title('Refuge Niche')
text(0.7,6.5,num2str(round(RfgNiche(1),1)))
text(1.8,11.5,num2str(round(RfgNiche(2),1)))
text(2.7,9.5,num2str(round(RfgNiche(3),1)))

% Optimum difference
subplot(2,3,6)
bar(DiffNiche,'FaceColor','k','EdgeColor','none','FaceAlpha',0.5)
box off
ylim([0 30])
xticklabels({'Mic-B/MicV-B','Mic-A/MicV-A','Mic-C/MicV-C'})
xtickangle(45)
title('Optimum difference')
text(0.7,5.5,num2str(round(DiffNiche(1),1)))
text(1.8,11.5,num2str(round(DiffNiche(2),1)))
text(2.7,8.5,num2str(round(DiffNiche(3),1)))

set(findall(gcf,'-property','FontSize'),'FontSize',16,'FontName','Times New Roman')

%% Print
print(hfig,'-depsc','-r600','outputs/Figure4_R0.eps')

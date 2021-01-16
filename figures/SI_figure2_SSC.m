%% A thermal trade-off between viral production and degradation drives phytoplankton-virus population dynamics
% === Code for: ===
% SI Figure 10 -- Side Scatter obtain by flow cytometry for Micromonas from
% Demory et al. ISME J. 2018
% Unpublished data © David Demory All Rights Reserved (contact me if needed)
% David Demory - Jan 2021

%% setup
data = xlsread('data/SSC.xlsx');
T = data(:,1);
SSC = data(:,2);

%% figure
hfig = figure('Color','white','position',[0 0 500 250])
plot(T,SSC,'k.','MarkerSize',15)
xlabel('Temperature ($^\circ$ C)','Interpreter','Latex')
ylabel('Average normalized SSC','Interpreter','Latex')

mdl = fitlm(T,SSC);
TT = 0:0.1:35;
m = 1.0376-0.0021131*TT;

hold on
a = plot(TT,m,'-','LineWidth',3)
c = get(a,'Color')
legend('data','model','orientation','horizontal','location','NorthWest')
legend boxoff
text(2,0.5,'y = -0.002 T + 1.04 -- R^2 = 0.01 -- p-value = 0.236','color',c)


% Save
print(hfig,'-depsc','-r600','outputs/SI_Figure2_SSC.eps')


%% A thermal trade-off between viral production and degradation drives phytoplankton-virus population dynamics
% === Code for: ===
% SI figure 11 -- Dynamics for varying temperature acording some scenarios
% David Demory - Jan 2021

%% Setup
% Cardinal parameters of MicB/MicV-B pair
Tmax = 33.46;
TR0  = 28.65;

% load hyperparameters
load('pmin_829_Arrhenius_v2.mat');

% create a figure
hfig = figure('color','White','position',[0 0 1000 500]);

% plots color for host and virus signals
colH = [26, 188, 156 ]/255;
colV  = [165, 105, 189]/255;

disp('Calculating the dynamics ... takes time!')

%% Case 1
% Average Temperature of the signal
T0 = TR0+0.2;
% Temperature amplitude = % of T0 
alpha = 0.1;
% Period of the signal 
f = 10;

% ---------- Case 1.1 = Temperature is fluctuating ----------
[t,vecT,tfit,yfit] = GenDyn(T0,alpha,f,pmin);

% Temperature subplot
han1 = subplot(2,4,1);
plot(t,vecT,'lineWidth',3);
yline(Tmax,'k:','LineWidth',3);
yline(TR0,'k-','LineWidth',3);
xlim([0 50]);
box off
ylabel('Temperature (^\circC)');
ylim([24 35]);
xlim([0 45]);
xticks([0,15,30,45]);
xticklabels({''});
% Dynamic subplot
subplot(2,4,5);
hold on
plot(tfit,yfit(:,1),'-','LineWidth',3,'Color',colH);
plot(tfit,yfit(:,2),':','LineWidth',3,'Color',colH);
plot(tfit,yfit(:,3),'-','LineWidth',3,'Color',colV);
plot(tfit,yfit(:,4),':','LineWidth',3,'Color',colV);
set(gca,'Yscale','log');
xlim([0 50]);
ylim([1 1E10]);
ylabel('Concentration (per ml)');
xlabel('Time (days)');
xlim([0 45]);
xticks([0,15,30,45]);

% ---------- Case 1.2 = Temperature is constant ----------
[t,vecT,tfit,yfit] = GenDyn(T0,0,f,pmin);

% Temperature subplot
han2 = subplot(2,4,2);
plot(t,vecT,'lineWidth',3);
yline(Tmax,'k:','LineWidth',3);
yline(TR0,'k-','LineWidth',3);
box off;
ylim([24 35]);
xlim([0 45]);
xticks([0,15,30,45]);
xticklabels({''});
% Dynamic subplot
subplot(2,4,6);
hold on
plot(tfit,yfit(:,1),'-','LineWidth',3,'Color',colH);
plot(tfit,yfit(:,2),':','LineWidth',3,'Color',colH);
plot(tfit,yfit(:,3),'-','LineWidth',3,'Color',colV);
plot(tfit,yfit(:,4),':','LineWidth',3,'Color',colV);
set(gca,'Yscale','log');
xlim([0 50]);
ylim([1 1E10]);
xlabel('Time (days)');
xlim([0 45]);
xticks([0,15,30,45]);

%% case 2
% Average Temperature of the signal
T0 = TR0-1;
% Temperature amplitude = % of T0 
alpha = 0.1;
% Period of the signal 
f = 10;

% ---------- Case 2.1 = Temperature is fluctuating ----------
[t,vecT,tfit,yfit] = GenDyn(T0,alpha,f,pmin);

% Temperature subplot
han1 = subplot(2,4,3);
plot(t,vecT,'lineWidth',3);
yline(Tmax,'k:','LineWidth',3);
yline(TR0,'k-','LineWidth',3);
xlim([0 50]);
box off
legend('Temperature','T_{max}','T_{R0 = 1}','location','SouthWestOutside','orientation','horizontal');
legend('boxoff');
ylim([24 35]);
xlim([0 45]);
xticks([0,15,30,45]);
xticklabels({''});
% Dynamic subplot
han3 = subplot(2,4,7);
hold on
plot(tfit,yfit(:,1),'-','LineWidth',3,'Color',colH);
plot(tfit,yfit(:,2),':','LineWidth',3,'Color',colH);
plot(tfit,yfit(:,3),'-','LineWidth',3,'Color',colV);
plot(tfit,yfit(:,4),':','LineWidth',3,'Color',colV);
set(gca,'Yscale','log');
xlim([0 50]);
ylim([1 1E10]);
legend('S','I','V_{i}','V_{ni}','location','NorthOutside','orientation','horizontal');
legend('boxoff');
xlabel('Time (days)');
xlim([0 45]);
xticks([0,15,30,45]);

% ---------- Case 2.2 = Temperature is constant ----------
[t,vecT,tfit,yfit] = GenDyn(T0,0,f,pmin);
% Temperature subplot
han2 = subplot(2,4,4);
plot(t,vecT,'lineWidth',3);
yline(Tmax,'k:','LineWidth',3);
yline(TR0,'k-','LineWidth',3);
xlim([0 50]);
box off
ylim([24 35]);
xlim([0 45]);
xticks([0,15,30,45]);
xticklabels({''});
% Dynamic subplot
han4 = subplot(2,4,8);
hold on
plot(tfit,yfit(:,1),'-','LineWidth',3,'Color',colH);
plot(tfit,yfit(:,2),':','LineWidth',3,'Color',colH);
plot(tfit,yfit(:,3),'-','LineWidth',3,'Color',colV);
plot(tfit,yfit(:,4),':','LineWidth',3,'Color',colV);
set(gca,'Yscale','log');
xlim([0 50]);
ylim([1 1E10]);
xlabel('Time (days)');
xlim([0 45]);
xticks([0,15,30,45]);

%% Clean plots

han1.Position = [han1.Position(1) han2.Position(2) han1.Position(3) han2.Position(4)];
han3.Position = [han3.Position(1) han4.Position(2) han3.Position(3) han4.Position(4)];

set(findall(gcf,'-property','FontSize'),'FontSize',14);

subplot(2,4,1);
text(-5,36,'a','FontSize',18);

subplot(2,4,2);
text(-5,36,'b','FontSize',18);

subplot(2,4,3);
text(-5,36,'c','FontSize',18);

subplot(2,4,4);
text(-5,36,'d','FontSize',18);

set(gcf,'Color','White');

%% Save
print(hfig,'-depsc','-r600','outputs/SI_Figure11_dynamics.eps');





%% A thermal trade-off between viral production and degradation drives phytoplankton-virus population dynamics
% === Code for: ===
% SI figure 8 -- Non-infectious virus (Vni) fraction against temperature
% David Demory - Jan 2021

%% setup

load('pmin_829_Arrhenius_v2.mat');

colH = [ 26, 188, 156 ]/255;
colV = [ 165, 105, 189]/255;
TR0 = 28.59;
YY = [];
T = [0:1:33];

%% ODE integration and calcul of Vni fraction as function of temperature
disp('Integrating the ODE for all the temperature ... takes time!')

for i = 1:length(T);
    
    Ti = T(i);
    
    t0 = 0:1:365;
    x0 = [1E3,0,1E4,0];
    pas = 0.1;
    t  = 1:pas:365;
    tT = 0:1:365;
    options = odeset('RelTol',1E-6,'AbsTol',1E-6,'Maxstep',0.1,'Events',@myEventsFcn);
    [tfit,yfit] = ode45(@SIV2_v2,t,x0,options,pmin,Ti);
    
    YY = [YY,yfit(end,4)./(yfit(end,4)+yfit(end,3))];
    
end

%% Figure
hfig = figure('color','white');
plot(T/TR0,YY,'LineWidth',3)
xline(1,'k:','LineWidth',3);
xlabel('T/T_{R0=1}')
ylabel('Fraction of V_{ni}')
box off

% fontsize
set(findall(gcf,'-property','FontSize'),'FontSize',14)

% Save
print(hfig,'-depsc','-r600','outputs/SI_Figure9_fractionVni.eps')

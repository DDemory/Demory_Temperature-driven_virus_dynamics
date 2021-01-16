function [t,vecT,tfit,yfit] = GenDyn(T0,alpha,f,theta)

% Function that integrete the model for a temperature signal that
% fluctuates with a sinusoidal shape.
% scenarii.

% -- Temperature amplitude:
% alpha = fractio of meanT apply to deltaT.

% T0 is the mean temperature of the temperature signal.
% theta are the hyperparameters 

%% Temperature signal
pas = 0.1;
t  = 0:pas:365;
x0 = [1E6,0,1E7,0];
vecT = T0-alpha*T0*sin((1/f)*2*pi*t);

%% Integration
options = odeset('RelTol',1E-6,'AbsTol',1E-6,'Maxstep',0.1,'Events',@myEventsFcn);
[tfit1,yfit1]   = ode45(@SIV2_v2_T,t,x0,options,theta,vecT,t);

yfit = yfit1;
tfit = tfit1;
yfittemp = yfit1;
tfittemp = tfit1;

%for i = 1:3
if tfit(end) < t(end)
    id = 1;
    while tfit(end)+pas<t(end)
        id = id +1;
        x00 = yfit(end,:);
        x00(x00 <=1) = 0;
        [tfittemp,yfittemp]   = ode45(@SIV2_v2_T,tfit(end)+pas:pas:t(end),x00,options,theta,vecT,t);
        yfit = [yfit;yfittemp];
        tfit = [tfit;tfittemp];
    end
end


end


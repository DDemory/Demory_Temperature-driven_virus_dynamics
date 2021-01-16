
function [position,isterminal,direction] = myEventsFcn(t,y,theta,vecT,timeT)
% Script to stop the ODE calculation when y<0 - to be used in the odeset
% options for ode solver using 'Events'
% David Demory - Jan 2021

position = [y(1)-1;y(2)-1;y(3)-1;y(4)-1]; % The value that we want to be zero
isterminal = [1,1,1,1];  % Halt integration 
direction = [-1;-1;-1;-1];   % The zero can be approached from either direction

end
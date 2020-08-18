function dy = SIV2_v2_H0(t,x,theta,T)
% ODE system for the basal model
% David Demory

% t     = time
% x     = variables
% theta = parameters
% T     = temperature

% Parameters
mu = theta(1);
K = theta(2);
phi = theta(3);
lambda = theta(4);
epsilon = theta(5);
beta = theta(6);
sigma = theta(7);
delta = theta(8);
omega = theta(9);
psi = theta(10);

% Variables
S = x(1);
I = x(2);
Vi = x(3);
Vni = x(4);
Vtot = Vi+Vni;

% ODE system
dy(1) = mu*S.*(1-(S+I)./K) - phi*S.*Vi -psi*S;
dy(2) = phi*S.*Vi - lambda*I- psi*I;
dy(3) = (1-epsilon)*beta*lambda*I - phi*S.*Vi - sigma*Vi - omega*Vtot.*Vi;
dy(4) = epsilon*beta*lambda*I + sigma*Vi - delta*Vni - omega*Vtot.*Vni;

dy = dy';
end

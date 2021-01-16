function dy = SIV2_v2_log(t,x,thetalog,T)
% ODE System of the temperature-driven model with log-parameters as input
% David Demory -- Jan 2021

theta=exp(thetalog);

% Calculate model paramters given the temperature-driven function hyper-parameters
[mu,K,phi,lambda,epsilon,beta,sigma,delta,omega,psi] = Tdriven_fct(T,theta);

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

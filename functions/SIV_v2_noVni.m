function dy = SIV_v2_noVni(t,x,theta,T)
% ODE System of the temperature-driven model
% David Demory

% t     = time
% x     = Variables
% theta = Hyper-arameters
% T     = temperature

% Calculate model paramters given the temperature-driven function hyper-parameters
[mu,K,phi,lambda,beta,delta,omega,psi] = Tdriven_fct_SIVnoVni(T,theta);

% Variables
S = x(1);
I = x(2);
Vi = x(3);

% ODE system
dy(1) = mu*S.*(1-(S+I)./K) - phi*S.*Vi -psi*S;
dy(2) = phi*S.*Vi - lambda*I- psi*I;
dy(3) = beta*lambda*I - phi*S.*Vi - delta*Vi - omega*Vi.^2;

dy = dy';
end

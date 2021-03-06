function [mu,K,phi,lambda,beta,delta,omega,psi] = Tdriven_fct(T,para)
% Temperature-driven function
% David Demory

% T = temperature
% para = hyper-parameter of the temperature-driven functions

%% --- Host growth ---
TK = T+273.15;
A1 = para(1);
E1 = para(2);
A2 = para(3);
E2 = para(4);
f = A1*exp(-E1./TK);
g = A2*exp(-E2./TK);
mu = f;
psi = g;

%% --- Carrying capacity ---
K = para(5);

%% --- Viral adsorption ---
phiK = para(6);
Tphi = para(7);
phir = para(8);
numA = phiK;
denA = phiK+exp(-phir.*(TK-Tphi));
phi = numA./denA;

%% --- Lysis rate ---
s1 = para(9);
d1 = para(10);
s2 = para(11);
d2 = para(12);
f = s1*exp(-d1./TK);
g = s2*exp(-d2./TK);
lambda = f - g;
lambda(lambda<0) = 0;

%% --- Burst size ---
b1 = para(13);
beta = b1*(f-g);
beta(beta<0) = 0;

%% --- viral degradation ---
%delta parameters
del1 = para(14);
del2 = para(15);
delta = del1*exp(-del2./TK);

%% --- viral aggregation ---
omega = para(16);

end

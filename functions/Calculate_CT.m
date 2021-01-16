function [Topt,Tmax,muopt,Tmin] = Calculate_CT(para)
% Estimation of the cardinal temperatures related to the Hinshelwood
% functions (From Grimmaud thesis 2016)
% David Demory - Jan 2021

%% Hinshelwood parameters
A1 = para(1);
E1 = para(2);
A2 = para(3);
E2 = para(4);

%% Conversion to cardinal parameters
% Optimal host temperature
Topt = (E1-E2)/log((A1*E1)/(A2*E2))-273.15;
% Maximal host temperature
Tmax = (E1-E2)/log(A1/A2)-273.15;
% Optimal host growth rate
theta = (E2-E1)/E1;
muopt = theta*A2*exp(-E2/(Topt+273.15));
% Minimal host temperature
epsilon = 0.05;
gamma = log(((E2-E1)*A2*epsilon)/(E1*A1));
Tmin = (Topt*E1/gamma)/(Topt-(E2/gamma));
end

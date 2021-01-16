function [Hdata,Htime,Vdata,Vtime,Hsd,Vsd,para,pmin,T] = load_Experiments_forCi(strain)
%Load data and best parameters for the estimation of the confidence
%intervals.
%strain = 829,451 or 834
% David Demory -- Jan 2021

Hdata = [];Htime = []; Vdata = []; Vtime = []; Hsd = []; Vsd = [];
if strain == 829
    
    T = [9.5 12.5 20 25 27.5 30];
    para.T = T;
    load('pmin_829_Arrhenius.mat');
    
    for i = 1:6;
        % host
        tempH = xlsread('data_B829_HOST.xlsx',i);
        Hdata = [Hdata,tempH(:,2)];
        Htime = [Htime,tempH(:,1)];
        Hsd = [Hsd,tempH(:,3)];
        % virus
        tempV = xlsread('data_B829_VIRUS.xlsx',i);
        Vdata = [Vdata,tempV(:,2)];
        Vtime = [Vtime,tempV(:,1)];
        Vsd = [Vsd,tempV(:,3)];    
    end
    
elseif strain == 451
    
    T = [12.5 20 25 27.5 30];
    para.T = T;
    load('pmin_451_Arrhenius.mat');
    
    for i = 1:5;
        % host
        tempH = xlsread('data_A451_HOST.xlsx',i);
        Hdata = [Hdata,tempH(:,2)];
        Htime = [Htime,tempH(:,1)];
        
        % virus
        tempV = xlsread('data_A451_VIRUS.xlsx',i);
        Vdata = [Vdata,tempV(:,2)];
        Vtime = [Vtime,tempV(:,1)];
    end
    
elseif strain == 834
    
    T = [12.5 20 25 27.5 30];
    para.T = T;
    load('pmin_834_Arrhenius.mat');
    
    for i = 1:5;
        % host
        tempH = xlsread('data_C834_HOST.xlsx',i);
        Hdata = [Hdata,tempH(:,2)];
        Htime = [Htime,tempH(:,1)];
        
        % virus
        tempV = xlsread('data_C834_VIRUS.xlsx',i);
        Vdata = [Vdata,tempV(:,2)];
        Vtime = [Vtime,tempV(:,1)];
  
    end
end


para.Htime = Htime;
para.Vtime = Vtime;
para.Hdata = Hdata;
para.Vdata = Vdata;

end

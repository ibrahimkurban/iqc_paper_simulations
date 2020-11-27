function param = find_params(L,m,method)
%%FIND_PARAMS 
% select parameters for the chosen method 
% which is given in Proposition 1 at pg 62 and Proposition 12 at pg 75
%
% param = find_params(L,m,method)
% param.alpha, param.beta, param.rho
% param.rho is theoretical convergence rate
%
kappa = L/m;
alpha = [];
beta = [];

% Proposition 1 at page 62
switch method
    
    %GD popular choice
    case 'gd_pop' 
        alpha = 1/L;
        rho = 1-1/kappa;
        
    %NAG standard choice
    case 'nag_std' 
        alpha = 1/L;
        beta = (sqrt(kappa)-1)/(sqrt(kappa)+1);
        rho = 1-1/sqrt(kappa);
        
    %GD optimal tuning
    case 'gd_opt' 
        alpha = 2/(L+m);
        rho = (kappa-1)/(kappa+1);
        
    %NAG optimal tuning
    case 'nag_opt'
        alpha = 4/(3*L+m);
        beta = (sqrt(3*kappa+1)-2)/(sqrt(3*kappa+1)+2);
        rho = 1-2/(sqrt(3*kappa+1));
        
    %HBM optimal tuning
    case 'hbm_opt' 
        alpha = 4/(sqrt(L)+sqrt(m))^2;
        beta = (sqrt(kappa)-1)/(sqrt(kappa)+1)^2;
        rho = (sqrt(kappa)-1)/(sqrt(kappa)+1);
        
end

param.alpha = alpha;
param.beta  = beta;
param.rho2   = rho;
end
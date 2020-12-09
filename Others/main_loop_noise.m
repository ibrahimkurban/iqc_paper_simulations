
clear all
close all
clc
%%
N       = 25; %number of kappa
tol     = 1e-2;
kappas  = logspace(1e-4, 4, N);
method  = 'nag_std';

%allocation
r       = zeros(N,6);
for i=1:N
    i
    % smoothnes and strong-cvx parameters
    L   = 1;
    m   = L/kappas(i);
    
    % Parameters of the method
    params = find_params(L,m,method);
    
    % system matrices of the chosen method
    G   = getG(params,method);
    
    % system matrices of the weigted off by one iqc
    psi = @(r2b)(weighted_iqc(m,L,r2b));
    
    % solve the lmi
    r(i,5)  = sqrt(bisection_cvx('noise', G,psi,tol, 0.5));
    r(i,4)  = sqrt(bisection_cvx('noise', G,psi,tol, 0.4));
    r(i,3)  = sqrt(bisection_cvx('noise', G,psi,tol, 0.3));
    r(i,2)  = sqrt(bisection_cvx('noise', G,psi,tol, 0.2));
    r(i,1)  = sqrt(bisection_cvx('noise', G,psi,tol, 0.1));
    r(i, 6) = params.rho;
end

%% figure
figure
semilogx(kappas, r', 'linewidth', 2)
grid on



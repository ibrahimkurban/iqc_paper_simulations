
clear all
close all
clc
%%
N       = 100; %number of kappa
tol     = 1e-2;
kappas  = logspace(1e-4, 5, N);
method  = 'nag_std';

%allocation
r2_1    = zeros(N,1);
r2_2    = zeros(N,1);
rho     = zeros(N,1);
for i=1:N
    
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
    r2_1(i)  = bisection_cvx('other', G,psi,tol);
    r2_2(i)  = bisection_cvx('naive', G,psi,tol);
    rho(i)   = params.rho;
end

%% figure
figure
semilogx(kappas, rho, 'linewidth', 2); hold on
semilogx(kappas, sqrt(r2_1), 'linewidth', 2)
semilogx(kappas, sqrt(r2_2), 'linewidth', 2)
grid on




clear all
close all
clc
%%
tol     = 1e-2;
kappas  = 1e2; logspace(1e-4, 5, 100);
method  = 'nag_std';

% smoothnes and strong-cvx parameters
L   = 1;
m   = L/kappas;

% Parameters of the method
params = find_params(L,m,method);

% system matrices of the chosen method
G   = getG(params,method);

% system matrices of the weigted off by one iqc
psi = @(r2b)(weighted_iqc(m,L,r2b));

% solve the lmi
r2_1  = bisection_cvx('noise', G,psi,tol, 0.4);
r2_2  = bisection_cvx('naive', G,psi,tol);


clear all
close all
clc
%%
delta = 0.05; % noise
tol     = 1e-2;
kappas  = 1e1; logspace(1e-4, 5, 100);
method  = 'gd_opt'; 'nag_std';

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
r_1  = sqrt(bisection_cvx('noise', G,psi,tol, delta));
r_2  = sqrt(bisection_cvx('noise naive', G,psi,tol, delta));

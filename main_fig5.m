% This script produces the Figure 5 at page 80

%% Prep
clear all
close all
clc

%% Interpreter
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

%% Method
N       = 31; %number of kappa
tol     = 1e-2;
kappas  = logspace(0, 3, N);
kappas(1) = 10^(1e-5);
method  = 'hbm_opt';

% allocation
r2_1    = zeros(N,1);
r2_2    = zeros(N,1);
rho     = zeros(N,1);
rho_opt = zeros(N,1);
rho_quad= zeros(N,1);
rho_gd  = zeros(N,1);

for i=1:N
    
    % smoothness and strong-cvx parameters
    L   = 1;
    m   = L/kappas(i);
    
    % Parameters of the method
    params = find_params(L,m,method);
    
    % system matrices of the chosen method
    G   = getG(params,method);
    
    % system matrices of the weigted off by one iqc
    psi = @(r2b)(weighted_iqc(m,L,r2b));
    
    % solve the LMI
    % Weighted off-by-one LMI conservative
    r2_1(i)  = bisection_cvx('naive', G,psi,tol);
    % Weighted off-by-one LMI less conservative
    r2_2(i)  = bisection_cvx('other', G,psi,tol);
    
    % Theoretical lower bound
    rho_opt(i)  = (sqrt(kappas(i))-1)/(sqrt(kappas(i))+1);
    % Optimal Gradient rate
    tmp         = find_params(L,m,'gd_opt');
    rho_gd(i)   = tmp.rho;
end

%% Plot the figure
figure;
% semilogx(kappas, sqrt(r2_1), 'linewidth', 2);
semilogx(kappas, sqrt(r2_2), 'linewidth', 2); hold on
semilogx(kappas, rho_gd, 'linewidth', 2);
semilogx(kappas, rho_opt, 'linewidth', 2);
semilogx(kappas, ones(size(kappas)), 'k--', 'linewidth', 1);
ylim([0 1.1])
grid on
xlabel('Condition ratio L/m')
ylabel('Convergence rate $\rho$')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
legend({'LMI (weighted off-by-one)',...
    'Optimal Gradient rate',...
    'Theoretical lower bound'})

%% Save results
save('fig5_results_v2.mat')
savefig('figure5_v2')

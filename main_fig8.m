% This script produces the Figure 8 at page 82

%% Prep
clear all
close all
clc

%% Interpreter
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

%% Method
N       = 21; %number of kappa
Nb      = 51; %number of betas for grid search
tol     = 1e-2;
kappas  = logspace(0, 4, N);
kappas(1) = 10^(1e-5);
betas   = linspace(0, 1, Nb);
method  = 'hbm_opt';

% allocation
r2_1_grid   = zeros(N,Nb);
r2_2_grid   = zeros(N,Nb);
r2_1        = zeros(N,1);
r2_2        = zeros(N,1);
rho         = zeros(N,1);
rho_opt     = zeros(N,1);
rho_quad    = zeros(N,1);
rho_gd_cvx  = zeros(N,1);
rho_gd_quad = zeros(N,1);

for i=1:N
    i
    % Smoothness and strong-cvx parameters
    L   = 1;
    m   = L/kappas(i);
    
    % system matrices of the weigted off by one iqc
    psi = @(r2b)(weighted_iqc(m,L,r2b));
        
    % Grid search for beta
    for j=1:Nb
        % Parameters of the method
        alpha = 1/L; % fixed alpha
        beta = betas(j); % varying beta
        
        % system matrix for HBM
        G   = genericG(alpha,beta,0);
    
        % solve the LMI
        % Weighted off-by-one LMI conservative
        r2_1_grid(i,j)  = bisection_cvx('naive', G,psi,tol);
        % Weighted off-by-one LMI less conservative
        r2_2_grid(i,j)  = bisection_cvx('other', G,psi,tol);
    end
    
    % Grid search optimal rates
    r2_1(i) = min(r2_1_grid(i,:));
    r2_2(i) = min(r2_2_grid(i,:));
    % Theoretical lower bound
    rho_opt(i)  = (sqrt(kappas(i))-1)/(sqrt(kappas(i))+1);
    % Gradient rate with alpha=1/L (popular choice)
    tmp         = find_params(L,m,'gd_pop');
    rho_gd_cvx(i)   = tmp.rho.c;
    rho_gd_quad(i)  = tmp.rho.q;
end

%% Plot the figure
figure;
semilogx(kappas, sqrt(r2_1), 'linewidth', 2); hold on
semilogx(kappas, sqrt(r2_2), 'linewidth', 2);
semilogx(kappas, rho_gd_cvx, 'linewidth', 2);
semilogx(kappas, rho_opt, 'linewidth', 2);
grid on
xlabel('Condition ratio L/m')
ylabel('Convergence rate $\rho$')
set(findall(gcf,'-property','FontSize'),'FontSize',14)
legend({'LMI (weighted off-by-one cons)',...
    'LMI (weighted off-by-one less cons)',...
    'Gradient rate with $\alpha=\frac{1}{L}$',...
    'Theoretical lower bound'})

%% Save results
save('fig8_results_v2.mat')
savefig('figure8_v2')

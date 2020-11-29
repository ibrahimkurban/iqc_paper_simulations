% This script produces the Figure 12 at page 85

%% Prep
clear all
close all
clc

%% Interpreter
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

%% Method
N       = 41; %number of kappa
tol     = 1e-2;
kappas  = logspace(0, 4, N);
kappas(1) = 10^(1e-5);
deltas  = [0.05,0.1,0.2,0.3,0.4,0.5];
Nd      = numel(deltas); %number of deltas
method  = 'nag_std';

% allocation
r       = zeros(N,Nd+1);

for i=1:N
    
    % Smoothness and strong-cvx parameters
    L   = 1;
    m   = L/kappas(i);
    
    % Parameters of the method
    params = find_params(L,m,method);
    
    % system matrices of the chosen method
    G   = getG(params,method);
    
    % system matrices of the weigted off by one iqc
    psi = @(r2b)(weighted_iqc(m,L,r2b));
    
    % solve the LMI for different values of delta
    for j=1:Nd
        r(i,j)  = sqrt(bisection_cvx('noise', G,psi,tol, deltas(j)));
    end
    % Standard Nesterov quadratics
    r(i,Nd+1)   = params.rho.q;
end

%% Plot the figure
legendstr = cell(1,Nd+1);
for j=1:Nd
    legendstr{j} = ['Standard Nesterov for $\delta$=' num2str(deltas(j))]; 
end
legendstr{Nd+1}  = 'Standard Nesterov (quadratics)';

figure;
semilogx(kappas, r', 'linewidth', 2); hold on;
semilogx(kappas, ones(size(kappas)), 'k--', 'linewidth', 1);
ylim([0 1.2])
grid on
xlabel('Condition ratio L/m')
ylabel('Convergence rate $\rho$')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
legend(legendstr)

%% Save results
save('fig12_results_v2.mat')
savefig('figure12_v2')
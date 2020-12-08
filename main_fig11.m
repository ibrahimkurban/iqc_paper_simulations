% This script produces the Figure 11 at page 84

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
deltas  = [0.01,0.02,0.05,0.1,0.2,0.5];
Nd      = numel(deltas); %number of deltas
method  = 'gd_pop';

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
    
    % Gradient rate with alpha=2/(L+m) (optimal)
    tmp         = find_params(L,m,'gd_opt');
    r(i,Nd+1)   = tmp.rho;
end

%% Plot the figure
legendstr = cell(1,Nd+1);
for j=1:Nd
    legendstr{j} = ['Rate for $\delta$=' num2str(deltas(j))]; 
end
legendstr{Nd+1}  = 'Noise-free Gradient rate';

figure;
semilogx(kappas, r(:,1:Nd)', 'linewidth', 2); hold on;
semilogx(kappas, r(:,Nd+1)', 'k-.', 'linewidth', 2);
semilogx(kappas, ones(size(kappas)), 'k--', 'linewidth', 1);
ylim([0 1.1])
grid on
xlabel('Condition ratio L/m')
ylabel('Convergence rate $\rho$')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
legend(legendstr)

%% Save results
save('fig11_results_v2.mat')
savefig('figure11_v2')
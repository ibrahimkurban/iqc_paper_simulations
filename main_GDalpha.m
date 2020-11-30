% This script tries to answer the following question at Page 76 
% "What range of stepsizes can yield a given rate?"

% Note:
% We used following Matlab package while plotting, please download from:
% https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar
% by Rob Campbell

%% Prep
clear all
close all
clc

%% Interpreter
set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

%% Method
kappas  = [1.2,5,10,20,50,100];
N       = numel(kappas); %number of kappa
rhos    = [0.1:0.1:0.7,0.75:0.05:0.9,0.91:1e-2:0.99,0.991:1e-3:0.999,0.9991:1e-4:0.9999];
Nr      = numel(rhos); %number of rhos

% allocation
alpha_min   = zeros(N,Nr);
alpha_max   = zeros(N,Nr);

for i=1:N

    % Smoothness and strong-cvx parameters
    L   = 1;
    m   = L/kappas(i);
        
    % GD LMI for each given rate
    for j=1:Nr
        
        % square of convergence rate
        r2 = rhos(j)^2;
        
        % solve the LMI
        alpha = GD_LMI_alpha(r2,L,m); 
            
        % min and max stepsizes
        alpha_min(i,j) = alpha.min;
        alpha_max(i,j) = alpha.max;
    end

end

%% Plot the figure
colors = {'-r','-b','-m','-c','-g','-y'};
legendstr = cell(N,1);
figure;
hold on;
for i=1:N
    legendstr{i} = ['$\kappa=' num2str(kappas(i)) '$'];
    ind = find(alpha_min(i,:)>0,1,'first');
    mean_alpha = (alpha_min(i,ind:Nr)+alpha_max(i,ind:Nr))/2;
    errbar = [alpha_max(i,ind:Nr)-alpha_min(i,ind:Nr);zeros(1,Nr-ind+1)];
    shadedErrorBar(rhos(ind:Nr),alpha_min(i,ind:Nr),errbar,'lineprops',colors{i},'transparent',1);
end
xlim([0 1])
grid on
xlabel('Convergence rate $\rho$')
ylabel('Step size $\alpha$')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
legend(legendstr)

%% Save results
save('GDalpha_results.mat')
savefig('figure_GDalpha')

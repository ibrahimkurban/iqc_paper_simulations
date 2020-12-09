function G = genericG(alpha,beta1,beta2)
%%GENERICG
% Parametrization of G with three parameters given at Page 88 Eq. 6.1 
% 
% Specific usage for three algorithms:
% (alpha,beta1,beta2) = (alpha,0,0) for the Gradient method
% (alpha,beta1,beta2) = (alpha,beta,0) for the Heavy-ball method
% (alpha,beta1,beta2) = (alpha,beta,beta) for the Nesterov's method
%
%   G = genericG(alpha,beta1,beta2)
%       where G.A, G.B, G.C. G.D

G.A = [beta1+1, -beta1; 1, 0];
G.B = [-alpha; 0];
G.C = [beta2+1, -beta2];
G.D = 0;

end
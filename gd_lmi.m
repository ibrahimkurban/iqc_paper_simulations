function LMI = gd_lmi(rho2, alpha, L, m, lam)
%%GD_LMI
% construct 3x3 lmi given at 76
% which is sector IQC applied on GD
% LMI is linear in all parameters
%   alpha: step size
%   L    : smoothness parameter
%   m    : strong convexity parameter
%   rho2 : square of worst case convergence rate
%   lam  : is a positive constant
%
%   LMI = gd_lmi(rho2, alpha, L, m, lam)
%
LMI = [-2*m*L*lambda - rho2, (L+m)*lam, 1;
       (L+m)*lam,            -2*lam,    -alpha;
       1,                    -alpha,    -1];
end
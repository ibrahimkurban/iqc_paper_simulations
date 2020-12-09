function alpha = GD_LMI_alpha(r2,L,m)
%%GD_LMI_ALPHA
% finds the range of stepsizes that yield a given convergence rate by using
% the LMI formulation for Gradient Descent given in Eq. 4.8 at Page 76
%
% alpha = GD_LMI_alpha(r2,L,m)
%   where alpha.min, alpha.max are the min and max stepsize given a rate
% r2 : square of convergence rate
% L  : smoothness parameter
% m  : strong convexity parameter

% minimum alpha
cvx_begin sdp quiet
    variable alph
    variable lam
    lmi = gd_lmi(r2,alph,L,m,lam);
    minimize alph
    subject to
    lmi <= 0;
cvx_end

switch cvx_status
    case 'Solved'
        alpha.min = alph;
    otherwise
        alpha.min = 0;
end

% maximum alpha
cvx_begin sdp quiet
    variable alph
    variable lam
    lmi = gd_lmi(r2,alph,L,m,lam);
    maximize alph
    subject to
    lmi <= 0;
cvx_end
switch cvx_status
    case 'Solved'
        alpha.max = alph;
    otherwise
        alpha.max = 0;
end

end

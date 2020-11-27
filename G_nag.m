function nag = G_nag(alpha,beta)
%%G_NAG
% construct dynamical system matrices for NAG
% which is given at pg 60 eq. 2.5
%
% nag = G_nag(alpha,beta)
%   where nag.A nag.B, nag.C nag.D
%
nag.A = [1+beta -beta; 1 0];
nag.B = [-alpha;0];
nag.C = [1+beta, -beta];
nag.D = 0;
end
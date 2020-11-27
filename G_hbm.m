function hbm = G_hbm(alpha,beta)
%%G_HBM
% construct dynamical system matrices for HBM
% which is given at pg 61 eq. 2.7
%
% hbm = G_hbm(alpha,beta)
%   where hbm.A hbm.B, hbm.C hbm.D
%
hbm.A = [1+beta -beta; 1 0];
hbm.B = [-alpha;0];
hbm.C = [1, 0];
hbm.D = 0;
end
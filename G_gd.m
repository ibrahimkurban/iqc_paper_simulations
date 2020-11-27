function gd = G_gd(alpha)
%%G_GD
% construct dynamical system matrices for GD
% which is given at pg 60 eq. 2.4
%
% gd = G_gd(alpha)
%   where gd.A gd.B, gd.C gd.D
%
gd.A = 1;
gd.B = -alpha;
gd.C = 1;
gd.D = 0;
end
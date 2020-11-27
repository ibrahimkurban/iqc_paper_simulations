function Gh = G_interconnect(G,psi)
%%G_INTERCONNECT
% construct system matrices for the interconnected system G and psi
% which is given at pg 67 eq 3.6
%
% Gh = G_interconnect(G,psi)
%   where input is psi.A, psi.Bu psi.By, psi.C. psi.Du psi.Dy 
%               G.A G.B, G.C G.D
%         output is Gh.A Gh.B, Gh.C Gh.D
%
n1 = size(G.A,1);
Gh.A = [G.A, zeros(n1,1);psi.By*G.C psi.A];
Gh.B = [G.B;psi.Bu];
Gh.C = [psi.Dy*G.C psi.C];
Gh.D = psi.Du;
end
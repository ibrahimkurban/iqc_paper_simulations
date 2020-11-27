function LMI = iqc_lmi_1(P,r2,lam, G,psi)
%%IQC_LMI_1
% construct LMI at pg 67
% formulation in 3.9 of Theorem 4
%
% LMI = iqc_lmi_1(P,r2,lam,G,psi)
%
Gh  = G_interconect(G, psi);

Z1  = zeros(size(Gh.B,1));
M   = [0 1;1 0];
LMI = [Gh.A Gh.B]'*P*[Gh.A Gh.B] - r2*blkdiag(P, Z1)...
        + lam*[Ch Dh]'*M*[Ch Dh];
end
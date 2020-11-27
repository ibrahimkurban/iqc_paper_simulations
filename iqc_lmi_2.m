function LMI = iqc_lmi_2(P,r2,lam1, lam2, G,psi)
%%IQC_LMI_2
% construct LMI at pg 73
% formulation in remark 11 of Lemma 10
%
% LMI = iqc_lmi_2(P,r2,lam1, lam2, G,psi)
%

Gh  = G_interconect(G, psi);

Z1  = zeros(size(Gh.B,1));
Z2  = zeros(size(Gh.B,1),1);
M   = [0 1;1 0];
LMI = [Ah Bh]'*P*[Ah Bh] - r2*blkdiag(P, Z1)...
      + lam2*[Ch Dh]'*M*[Ch Dh]...
      + lam1*[psi.Dy*G.C, Z2, Dh]'*M*[psi.Dy*G.C, Z2, Dh];
end
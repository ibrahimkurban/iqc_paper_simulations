function LMI = iqc_lmi_2(P,r2,lam1, lam2, G,psi)
%%IQC_LMI_2
% construct LMI at pg 73
% formulation in remark 11 of Lemma 10
%
% LMI = iqc_lmi_2(P,r2,lam1, lam2, G,psi)
%

Gh  = G_interconnect(G, psi);

Z1  = zeros(size(Gh.B,2));
Z2  = zeros(size(Gh.D,1),1);
M   = [0 1;1 0];
LMI = [Gh.A Gh.B]'*P*[Gh.A Gh.B] - r2*blkdiag(P, Z1)...
      + lam2*[Gh.C Gh.D]'*M*[Gh.C Gh.D]...
      + lam1*[psi.Dy*G.C, Z2, Gh.D]'*M*[psi.Dy*G.C, Z2, Gh.D];
end
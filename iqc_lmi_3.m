function LMI = iqc_lmi_3(P,r2,lam0, lam1, lam2, G,psi,delta)
%%IQC_LMI_3
% construct LMI at pg 81
% formulation in eq 5.1 and use the trick in remark 11 at pg 73
%
% LMI = iqc_lmi_3(P,r2,lam0, lam1, lam2, G,psi,delta)
%
Gh  = G_interconnect_noise(G, psi);


Z1  = zeros(size(Gh.B,2));
Z2  = zeros(size(Gh.D,1),1);
Z3  = zeros(2,size(Gh.A,1));
M1  = [0 1;1 0]; 
M2  = [delta^2-1 1;1 -1];
LMI = [Gh.A Gh.B]'*P*[Gh.A Gh.B] - r2*blkdiag(P, Z1)...
        + lam2*[Gh.C Gh.D]'*M1*[Gh.C Gh.D]...
        + lam1*[psi.Dy*G.C, Z2, Gh.D]'*M1*[psi.Dy*G.C, Z2, Gh.D]...
        + lam0*[Z3, M1]'*M2*[Z3, M1];
end
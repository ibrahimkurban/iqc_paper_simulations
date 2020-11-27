function LMI = iqc_lmi_4(P,r2,lam1, lam2, G,psi,delta)
%%IQC_LMI_4
% construct LMI at pg 81
% formulation in eq 5.1
%
% LMI = iqc_lmi_4(P,r2,lam1, lam2, G,psi,delta)
%
Gh  = G_interconnect_noise(G, psi);


Z1  = zeros(size(Gh.B,2));
Z2  = zeros(2,size(Gh.A,1));
M1  = [0 1;1 0]; 
M2  = [delta^2-1 1;1 -1];
LMI = [Gh.A Gh.B]'*P*[Gh.A Gh.B] - r2*blkdiag(P, Z1)...
        + lam1*[Gh.C Gh.D]'*M1*[Gh.C Gh.D]...
        + lam2*[Z2, M1]'*M2*[Z2, M1];
end
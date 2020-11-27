function LMI = iqc_lmi_3(P,r2,lam1, lam2, G,psi,delta)
%%IQC_LMI_3
% construct LMI at pg 81
% formulation in eq 5.1
%
% LMI = iqc_lmi_3(P,r2,lam1, lam2, G,psi,delta)
%
Gh  = G_interconnect(G, psi);

Z1  = zeros(size(Gh.B,1),1);
Z2  = zeros(size(Gh.B,2));
Z3  = zeros(size(Gh.D,1),1);
M1  = [0 1;1 0]; 
M2  = [delta^2-1 1;1 -1];
LMI = [Gh.A Z1 Gh.B]'*P*[Gh.A Z1 Gh.B] - r2*blkdiag(P, Z2, Z2)...
        + lam1*[Gh.C Z3 Gh.D]'*M1*[Gh.C Z3 Gh.D]...
        + lam2*[0 0 0 0 1;0 0 0 1 0]'*M2*[0 0 0 0 1;0 0 0 1 0];
end
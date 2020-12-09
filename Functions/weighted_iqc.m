function psi = weighted_iqc(m,L,r2)
%%WEIGHTED_IQC
% construct Phi matrix by using weighted off by one IQC
% Lemma 10 at pg 72
%
% psi = weighted_iqc(m,L,r2)
%       where psi.A, psi.Bu psi.By, psi.C. psi.Du psi.Dy
%

psi.A   = 0;
psi.By  = -L;
psi.Bu  = 1;
psi.C   = [r2;0];
psi.Dy  = [L;-m];
psi.Du  = [-1; 1];
end
function [r2, P] = bisection_cvx(iqc_type, G,psi,tol,delta)
%%BISECTION_CVX
% solves the LMI given in
% iqc_type: 'naive'         formulation in 3.9 of Theorem 4
%            anything       formulation in remark 11 of Lemma 10 pg 73
% by using bisection method
%
% [r2, P] = bisection_cvx(iqc_type, G,psi,tol)
% r2 : square of convergence rate
%
%
%size of the problem
n    = size(G_interconnect(G,psi(1)).A,1);

%bisection interval
r2_h = 2;
r2_l = 0;
while(r2_h - r2_l>tol)
    %rho update
    r2 = 0.5*(r2_h + r2_l);
    %% cvx
    switch iqc_type
        case 'naive'
            lam = 1; % homogenity
            cvx_begin sdp quiet
            variable P(n,n) semidefinite
            lmi = iqc_lmi_1(P,r2,lam, G,psi(r2));
            minimize 0
            subject to
            lmi <= 0;
            cvx_end
        case 'noise'
            cvx_begin sdp quiet
            variable P(n,n) semidefinite
            variable lam1 nonnegative
            variable lam2 nonnegative
            lmi = iqc_lmi_3(P,r2,lam1, lam2, G,psi(r2), delta);
            minimize 0
            subject to
            lmi <= 0;
            lam1 + lam2 == 1; %homogenity
            cvx_end
        otherwise
            cvx_begin sdp quiet
            variable P(n,n) semidefinite
            variable lam1 nonnegative
            variable lam2 nonnegative
            lmi = iqc_lmi_2(P,r2,lam1, lam2,G,psi(1));
            minimize 0
            subject to
            lmi <= 0;
            lam1 + lam2 == 1; %homogenity
            lam2 <= r2*(lam1+lam2);
            cvx_end
    end
    %% bisection interval update
    switch cvx_status
        case 'Solved'
            r2_h = r2;
        otherwise
            r2_l = r2;
    end
end

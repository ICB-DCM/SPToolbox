% col_w_diff.m computes the Cholesky factorisation of a positive 
% definite matrix, A = L*L'. Furthermore, the derivatives of L are
% determined based on the derivatives of A.
%
% This code has been adapted from the paper:
%   Mike Giles, An extended collection of matrix derivative results for
%   forward and reverse mode algorithmic differentiation, Report no. 08/01 
%   from the Oxford University Computing Laboratory, Parks Road, Oxford, U.K.
%   url: http://people.maths.ox.ac.uk/gilesm/files/NA-08-01.pdf

function [L,dL] = chol_w_diff(A,dA)

N = size(A,1);
q = size(dA,3);

L = zeros(N,N);
dL = zeros(N,N,q);

for m = 1:N
    for n = 1:m
        for k = 1:n-1
            A(m,n) = A(m,n) - L(m,k)*L(n,k);
            for l = 1:q
                dA(m,n,l) = dA(m,n,l) - dL(m,k,l)*L(n,k) - L(m,k)*dL(n,k,l);
            end
        end
        if m==n
            L(m,m) = sqrt(A(m,m));
            for l = 1:q
                dL(m,m,l) = 0.5*dA(m,m,l)/L(m,m);
            end
        else
            L(m,n) = A(m,n)/L(n,n);
            for l = 1:q
                dL(m,n,l) = (dA(m,n,l)-L(m,n)*dL(n,n,l))/L(n,n);
            end
        end
    end
end
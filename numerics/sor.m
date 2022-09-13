% Solve matrix system using successive over-relaxation
%   [C] = sor(A,C,b,omega,tol)
%           A*C = b
% C provided must be initial guess.

function [C] = sor(A,C,b,omega,tol,max_iter)

    if ~exist('max_iter','var'), max_iter = 1000; end
   
    %First normalize the matrix such that all the diagonals are equal to 1. 
    for i=1:length(A)
        b(i)=b(i)/A(i,i);
        diag=A(i,i);
        for j=1:length(A)
            A(i,j)=A(i,j)/diag;
        end
    end
    k = 1;
    while abs(mean(A*C-b)) > tol && k <= max_iter
        for i = 1:length(C)
            sigma = 0;
            for j=1:length(C)
                if j ~= i
                    sigma = sigma + A(i,j)*C(j);
                end
            end
            C(i) = (1-omega)*C(i) + omega/A(i,i) * (b(i) - sigma);
        end
        k = k+1;
    end
    
    k
    if k >= max_iter, warning('Maximum iterations reached. Giving up.'); end
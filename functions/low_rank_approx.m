function [A_k, k] = low_rank_approx(A, tol)
%LOW_RANK_APPROX  Low-rank approximation via Singular Value Decomposition.
%
%   [A_k, k] = LOW_RANK_APPROX(A, tol)
%
%   INPUT:
%       A   - matrix to be approximated
%       tol - tolerance on the relative error in the Frobenius norm
%
%   OUTPUT:
%       A_k - rank-k approximation of the matrix A
%       k   - minimum rank such that the relative error is < tol
    [U, S, V] = svd(A);
    s = diag(S);
    n = length(s);
    
    % Frobenius norm of the original matrix (sqrt of the sum of squared singular values)
    normaA_F = sqrt(sum(s.^2));
    
    k = n; % Default value if a lower k is not found
    for i = 1:n-1
        % Error norm (singular values discarded from rank i+1 onwards)
        err_frob = sqrt(sum(s(i+1:end).^2));
        
        % Relative error
        if (err_frob / normaA_F) < tol
            k = i;
            break;
        end
    end
    
    fprintf("The optimal rank for the tolerance %e is: %d\n", tol, k);
    
    % Reconstruction of the approximated matrix Ak
    A_k = U(:, 1:k) * S(1:k, 1:k) * V(:, 1:k)';
end
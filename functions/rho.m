function [rhoGS, rhoJ] = rho(A)
% RHO  Calculates the spectral radius of the Jacobi and Gauss-Seidel iteration matrices.
%
%   [rhoGS, rhoJ] = RHO(A)
%
%   INPUT:
%       A     - square matrix of the linear system Ax = b
%
%   OUTPUT:
%       rhoGS - spectral radius of the Gauss-Seidel iteration matrix
%       rhoJ  - spectral radius of the Jacobi iteration matrix
%
%   The function builds the iteration matrices:
%       BJ  = I - D\A        (Jacobi)
%       BGS = I - T\A        (Gauss-Seidel)
%   where D is the diagonal part of A and T is the lower triangular part.
%
%   The spectral radius is calculated as:
%       rho = max(abs(eig(B)))
%
%   Based on the value of rho (< 1 or >= 1), the function prints
%   a message indicating whether the corresponding method converges.
    n = size(A, 1);
    I = eye(n);
    
    % Jacobi calculation
    D = diag(diag(A));
    BJ = I - D\A;
    rhoJ = max(abs(eig(BJ)));
    if rhoJ < 1
        fprintf("The Jacobi method converges with spectral radius %f\n", rhoJ);
    else
        fprintf("The Jacobi method does not converge with spectral radius %f\n", rhoJ);
    end
    
    % Gauss-Seidel calculation
    T = tril(A);
    BGS = I - T\A;
    rhoGS = max(abs(eig(BGS)));
    if rhoGS < 1
        fprintf("The Gauss-Seidel method converges with spectral radius %f\n", rhoGS);
    else
        fprintf("The Gauss-Seidel method does not converge with spectral radius %f\n", rhoGS);
    end
end
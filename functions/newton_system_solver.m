function [xvect, iter] = newton_system_solver(x0, tol, nmax, F, J)
% NEWTON_SYSTEM_SOLVER  Solves a system of non-linear equations using Newton's method.
%
%   [xvect, iter] = newton_system_solver(x0, tol, nmax, F, J)
%
%   INPUT:
%       x0   - initial vector (initial estimate of the solution)
%       tol  - tolerance on the step size (stopping criterion)
%       nmax - maximum number of iterations
%       F    - function handle that evaluates the system F(x)
%       J    - function handle that evaluates the Jacobian J(x)
%
%   OUTPUT:
%       xvect - matrix containing all iterate vectors column by column
%       iter  - number of iterations performed
err = tol + 1;
x = x0(:); 
xvect = x;
iter = 0;
while (err > tol && iter < nmax)
    F_val = F(x);
    J_val = J(x);
    delta = J_val\(-F_val);
    x = x + delta;
    xvect = [xvect, x];
    iter = iter + 1;
    err = norm(delta);
end
end
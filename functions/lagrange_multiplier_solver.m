function [lambda_sol, z, iter] = lagrange_multiplier_solver(phi, g, n, z0, tol, max_it)
% LAGRANGE_MULTIPLIER_SOLVER  Solves constrained optimization problems using 
%                   Lagrange multipliers and Newton's method.
%
%   [lambda_sol, z, iter] = lagrange_multiplier_solver(phi, g, n, z0, tol, max_it)
%
%   INPUT:
%       phi     - symbolic objective function to be optimized
%       g       - vector of symbolic constraints g_i(x) = 0
%       n       - number of decision variables
%       z0      - initial guess vector for Newton's method
%                 (contains both variables x and multipliers lambda)
%       tol     - tolerance on the norm of the increment
%       max_it  - maximum number of allowed iterations
%
%   OUTPUT:
%       lambda_sol - solution values for the Lagrange multipliers
%       z          - approximate solution for the decision variables
%       iter       - actual number of iterations performed
%      
%
%   DESCRIPTION:
%       The function constructs the Lagrangian:
%
%           L(x, lambda) = phi(x) + sum_i lambda_i * g_i(x)
%
%       and symbolically generates:
%           - the system F = ∇L = 0 (stationarity conditions)
%           - the Jacobian matrix J = ∂F/∂(x,lambda)
%
%       To ensure efficiency during the iterative loop, F and J are
%       converted into numerical functions using matlabFunction.
%
%       Newton's method is then applied:
%
%           J(z_k) * delta = -F(z_k)
%           z_{k+1} = z_k + delta
%
%       The loop terminates when:
%           norm(delta, inf) < tol
%       or when max_it is reached.
%
%       This function provides a general approach to solve optimization 
%       problems with any number of variables and constraints, combining 
%       symbolic flexibility with numerical efficiency.

% 1. Symbolic setup (executed once)
x = sym('x',[n,1]);
m = length(g);
lambda = sym('lambda',[m,1]);

L = simplify(phi + lambda.' * g); % Lagrangian
vars = [x; lambda];

F_sym = gradient(L, vars);       % System of equations
J_sym = jacobian(F_sym, vars);   % Jacobian matrix

% Conversion to numerical functions for speed inside the loop
F_num = matlabFunction(F_sym, 'Vars', {vars});
J_num = matlabFunction(J_sym, 'Vars', {vars});

% 2. Newton's method initialization
z = z0(:);
iter = 0;
err = tol + 1; % Initial value to enter the loop

% 3. Newton's loop with WHILE
while (err > tol) && (iter < max_it)
    F_val = F_num(z);
    J_val = J_num(z);

    % Solving the linear system J * delta = -F
    delta = J_val \ (-F_val);

    % Solution update
    z = z + delta;

    % Error computation (norm of the increment)
    err = norm(delta, inf); 
    iter = iter + 1;
end

lambda_sol = z(n+1:end);
z = z(1:n);

% Final check on convergence
if iter == max_it && err > tol
    warning('Maximum number of iterations reached without convergence.');
end
end
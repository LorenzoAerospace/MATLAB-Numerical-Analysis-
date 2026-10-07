function [U_com, x] = mixed_poisson_1d(f_fun, L, h, alpha, g_beta, mu, eta, sigma_fun)
% MIXED_POISSON_1D  Solves the variable-coefficient Advection-Diffusion-Reaction problem:
%                   -mu*u''(x) + eta*u'(x) + sigma(x)*u(x) = f(x) in (0, L)
%                   with mixed Dirichlet (left) and Neumann (right) boundary conditions.
%
%   [U_com, x] = MIXED_POISSON_1D(f_fun, L, h, alpha, g_beta, mu, eta, sigma_fun)
%
%   INPUT:
%       f_fun     - anonymous function for the source term, e.g., @(x) sin(x)
%       L         - length of the spatial domain [0, L]
%       h         - spatial discretization step size
%       alpha     - Dirichlet boundary condition at the left boundary: u(0) = alpha
%       g_beta    - Neumann boundary condition value at the right boundary: u'(L) = g_beta
%       mu        - diffusion coefficient (positive constant)
%       eta       - advection/velocity coefficient (signed constant)
%       sigma_fun - anonymous function for the reaction coefficient, e.g., @(x) x.^2
%
%   OUTPUT:
%       U_com     - column vector containing the complete solution (internal nodes + boundaries)
%       x         - column vector of the spatial grid nodes

% 1. Definition of the spatial grid
N_tot = round(L / h);     
x = (0:N_tot)' * h;     
N_unknowns = N_tot;    
x_inc = x(2:end); % Nodes associated with the unknowns (from x_1 to x_N_tot)

% 2. Initialization of the matrix and the right-hand side vector
% NOTE: sigma_fun(x_inc) returns a vector of the same length as x_inc.
% We use element-wise operations (.* and + between vectors)
main  = 2 * mu * ones(N_unknowns, 1) + sigma_fun(x_inc) * h^2; 
sopra = (-mu + eta * h / 2) * ones(N_unknowns - 1, 1);
sotto = (-mu - eta * h / 2) * ones(N_unknowns - 1, 1);
rhs = f_fun(x_inc) * h^2; 

% 3. Treatment of the Left Boundary (Dirichlet: u_0 = alpha)
coeff_u0 = -mu - (eta * h / 2);
rhs(1) = rhs(1) - coeff_u0 * alpha; 

% 4. Treatment of the Right Boundary (Centered Neumann with a ghost node)
% Modification of the sub-diagonal coefficient of the last row:
sotto(end) = -2 * mu; 
% Modification of the last element of the right-hand side vector:
rhs(end) = rhs(end) - 2 * h * g_beta * (-mu + (eta * h / 2));

% 5. Assembly of the tridiagonal matrix A
A = diag(main) + diag(sopra, 1) + diag(sotto, -1);

% 6. Resolution of the linear system
U_inc = A \ rhs;

% 7. Reconstruction of the complete solution vector
U_com = [alpha; U_inc];
end
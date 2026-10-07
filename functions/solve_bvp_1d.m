function [x, U] = solve_bvp_1d(mu, eta, sigma_fun, f_fun, alpha, gamma, L, right_bc_type, h)
% SOLVE_BVP_1D  Solves -mu*u'' + eta*u' + sigma(x)*u = f(x) on the interval [0, L]
%
%   [x, U] = SOLVE_BVP_1D(mu, eta, sigma_fun, f_fun, alpha, gamma, L, right_bc_type, h)
%
%   INPUT:
%       mu            - diffusion coefficient (positive constant)
%       eta           - advection/velocity coefficient (signed constant)
%       sigma_fun     - function handle for the reaction term, e.g., @(x) x.^2
%       f_fun         - function handle or scalar for the source term
%       alpha         - Dirichlet boundary condition value at the left boundary: u(0) = alpha
%       gamma         - boundary condition value at the right boundary (x = L)
%       L             - right endpoint of the spatial domain [0, L]
%       right_bc_type - string specifying the right boundary type:
%                       * 'Dirichlet'        - u(L) = gamma
%                       * 'Neumann_indietro' - u'(L) = gamma (Backward approximation)
%                       * 'Neumann_centrato' - u'(L) = gamma (Centered approximation with ghost node)
%       h             - spatial discretization step size
%
%   OUTPUT:
%       x             - column vector of spatial grid nodes
%       U             - column vector containing the complete solution (nodes 0 to M)

% 1. Grid construction
x = (0:h:L)';
M = length(x) - 1; % Number of intervals (Dimension of the internal system)

% Evaluation of functions in the nodes involved in the system (from x_1 to x_M)
% Note: x(1) is x_0=0, x(end) is x_M=L
s = sigma_fun(x);

% 2. Assembly of diagonals for internal nodes
main = (2 * mu / h^2 + s(2:end)); 
up   = (-mu / h^2 + eta / (2 * h)) * ones(M, 1);
down = (-mu / h^2 - eta / (2 * h)) * ones(M, 1);

A = sparse(1:M, 1:M, main, M, M) + ...
    sparse(2:M, 1:M-1, down(1:M-1), M, M) + ...
    sparse(1:M-1, 2:M, up(1:M-1), M, M);

% Right-hand side initialization
% If f is a function handle, evaluate it; if it is a constant scalar, create a vector
if isa(f_fun, 'function_handle')
    rhs = f_fun(x(2:end));
else
    rhs = f_fun * ones(M, 1);
end

% 3. Treatment of the left boundary (Dirichlet on u_0)
% Modify the first element of rhs with the 'down' coefficient multiplying alpha
val_down_bordo = -mu / h^2 - eta / (2 * h);
rhs(1) = rhs(1) - val_down_bordo * alpha;

% 4. Treatment of the right boundary (x = L, row M of the matrix)
switch lower(right_bc_type)
    case 'dirichlet'
        A(end, :) = 0;
        A(end, end) = 1;
        rhs(end) = gamma;

    case {'neumann_indietro', 'neumann_backward'}
        A(end, :) = 0;
        A(end, end) = 1 / h;
        A(end, end-1) = -1 / h;
        rhs(end) = gamma;

    case {'neumann_centrato', 'neumann_centered'}
        A(end, :) = 0;
        A(end, end-1) = -2 * mu / h^2;
        A(end, end)   = 2 * mu / h^2 + s(end);
        rhs(end) = rhs(end) - (-2 * mu / h + eta) * gamma;

    otherwise
        error('Unrecognized right boundary type!');
end

% 5. Resolution of the linear system and inclusion of u_0 = alpha
U_sol = A \ rhs;
U = [alpha; U_sol];
end
function [U_com, X, Y] = poisson_2d_dirichlet(f_fun, L_x, L_y, h, g_bordo, mu, eta_x, eta_y, sigma_fun)
% POISSON_2D_DIRICHLET  Solves the 2D problem with spatially variable reaction:
%                      -mu*Delta(u) + eta_x*u_x + eta_y*u_y + sigma(x,y)*u = f(x,y)
%                      with pure Dirichlet boundary conditions.
%
%   [U_com, X, Y] = POISSON_2D_DIRICHLET(f_fun, L_x, L_y, h, g_bordo, mu, eta_x, eta_y, sigma_fun)
%
%   INPUT:
%       f_fun     - anonymous function for the source term, e.g., @(x,y) sin(x.*y)
%       L_x       - length of the spatial domain along the x-axis [0, L_x]
%       L_y       - length of the spatial domain along the y-axis [0, L_y]
%       h         - spatial discretization step size
%       g_bordo   - Dirichlet boundary condition value (uniform on the boundary)
%       mu        - diffusion coefficient (positive constant)
%       eta_x     - advection/velocity coefficient along the x-axis
%       eta_y     - advection/velocity coefficient along the y-axis
%       sigma_fun - anonymous function for the reaction coefficient, e.g., @(x,y) x.^2 + y.^2
%
%   OUTPUT:
%       U_com     - complete solution matrix (Ny_tot+1 by Nx_tot+1) including boundaries
%       X         - 2D grid coordinates matrix along the x-axis
%       Y         - 2D grid coordinates matrix along the y-axis

% 1. Grid construction
Nx_tot = round(L_x / h); 
Ny_tot = round(L_y / h);
Nx = Nx_tot - 1;       
Ny = Ny_tot - 1; 
x = (0:Nx_tot)' * h;     
y = (0:Ny_tot)' * h;
[X, Y] = meshgrid(x, y); 

% Internal nodes
X_int = X(2:end-1, 2:end-1);
Y_int = Y(2:end-1, 2:end-1);

% 2. Construction of base 1D operators
Ix = speye(Nx); 
Iy = speye(Ny);

% Operators along X
t_x = zeros(Nx, 1); t_x(1) = 2; t_x(2) = -1;
Tx_diff = toeplitz(t_x);
cx = zeros(Nx, 1); if Nx > 1, cx(2) = 0.5; end
Tx_trasp = toeplitz(cx, -cx);

% Operators along Y
t_y = zeros(Ny, 1); t_y(1) = 2; t_y(2) = -1;
Ty_diff = toeplitz(t_y);
cy = zeros(Ny, 1); if Ny > 1, cy(2) = 0.5; end
Ty_trasp = toeplitz(cy, -cy);

% 3. Assembly of Diffusion and Advection contributions
A_diff    = mu * (kron(Iy, Tx_diff) + kron(Ty_diff, Ix));          
A_trasp_x = eta_x * h * kron(Iy, Tx_trasp);                       
A_trasp_y = eta_y * h * kron(Ty_trasp, Ix);                       

% --- MODIFICATION FOR VARIABLE SIGMA(X,Y) ---
% Evaluate sigma_fun on internal nodes (returns a Ny x Nx matrix)
sigma_val = sigma_fun(X_int, Y_int);
% Convert to a column vector to linearize nodes in the same order as f_val(:)
sigma_vett = sigma_val(:);
% Create a sparse diagonal matrix with the point-by-point sigma values
A_reaz = h^2 * spdiags(sigma_vett, 0, Nx*Ny, Nx*Ny);

% Total system matrix
A = A_diff + A_trasp_x + A_trasp_y + A_reaz;

% 4. Construction of the right-hand side vector (RHS)
f_val = f_fun(X_int, Y_int);
rhs = f_val(:) * h^2; 

% 5. Dirichlet boundary correction (if g_bordo ~= 0)
if g_bordo ~= 0
    % The reaction acts only on internal nodes and does not affect boundary exchange.
    % We use only the diffusive and advective parts to detect boundaries.
    A_bordo_detect = A_diff + A_trasp_x + A_trasp_y;
    Bordi = -sum(A_bordo_detect, 2); 
    rhs = rhs + Bordi * g_bordo;
end

% 6. Resolution of the sparse linear system
U_interni_vett = A \ rhs;

% 7. Reconstruction of the solution matrix (Ny x Nx)
U_interni = reshape(U_interni_vett, Ny, Nx);

% 8. Construction of the complete solution matrix including known boundaries
U_com = g_bordo * ones(Ny_tot+1, Nx_tot+1);
U_com(2:end-1, 2:end-1) = U_interni;
end
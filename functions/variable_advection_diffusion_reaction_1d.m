function [U_com, x] = variable_advection_diffusion_reaction_1d(mu, eta_fun, sigma_fun, f_fun, L, h, alpha, beta)
% VARIABLE_ADVECTION_DIFFUSION_REACTION_1D  Solves -mu*u'' + eta(x)*u' + sigma(x)*u = f(x)
%                                          with Dirichlet boundary conditions u(0) = alpha, u(L) = beta.
%
%   [U_com, x] = VARIABLE_ADVECTION_DIFFUSION_REACTION_1D(mu, eta_fun, sigma_fun, f_fun, L, h, alpha, beta)
%
%   INPUT:
%       mu        - diffusion coefficient (constant > 0)
%       eta_fun   - anonymous function for advection, e.g., @(x) 6 + 0*x or @(x) x
%       sigma_fun - anonymous function for reaction, e.g., @(x) 2*(4-x)
%       f_fun     - anonymous function for the forcing term, e.g., @(x) 12 + 0*x
%       L         - length of the spatial domain [0, L]
%       h         - spatial discretization step size
%       alpha     - Dirichlet condition at the left boundary u(0) = alpha
%       beta      - Dirichlet condition at the right boundary u(L) = beta
%
%   OUTPUT:
%       U_com     - column vector containing the complete solution (internal nodes + boundaries)
%       x         - column vector of spatial grid nodes

    % 1. Calculation of grid parameters
    N_tot = round(L / h);     % Total number of intervals
    N = N_tot - 1;          % Number of internal nodes
    x = (0:N_tot)' * h;     % Vector of all nodes from x_0 to x_N_tot
    x_int = x(2:end-1);     % Vector of internal nodes only (from x_1 to x_N)
    
    % Evaluation of coefficients at the internal nodes
    eta_val = eta_fun(x_int);
    sigma_val = sigma_fun(x_int);
    f_val = f_fun(x_int);
    
    % 2. Initialization of the diagonals of matrix A
    main  = zeros(N, 1);
    sopra = zeros(N-1, 1);
    sotto = zeros(N-1, 1);
    
    % 3. Row-by-row assembly (multiplying everything by h^2 to scale the RHS)
    % Common contribution of Diffusion (-mu * u''):
    % sub-diag: -mu, diagonal: 2*mu, super-diag: -mu
    % Loop over internal nodes to handle variable advection and adaptive Upwind
    for j = 1:N
        % Local reaction degree
        reaz_local = sigma_val(j) * h^2;
        
        % Determination of the local advection scheme (Upwind)
        if eta_val(j) >= 0
            % Wind towards the RIGHT (eta > 0) -> Backward Upwind
            % u' approx (u_j - u_{j-1})/h -> multiplied by h^2 becomes eta*h*(u_j - u_{j-1})
            trasp_main  = eta_val(j) * h;
            trasp_sotto = -eta_val(j) * h;
            trasp_sopra = 0;
        else
            % Wind towards the LEFT (eta < 0) -> Forward Upwind
            % u' approx (u_{j+1} - u_j)/h -> multiplied by h^2 becomes eta*h*(u_{j+1} - u_j)
            trasp_main  = -eta_val(j) * h;
            trasp_sotto = 0;
            trasp_sopra = eta_val(j) * h;
        end
        
        % Construction of the element on the main diagonal
        main(j) = 2 * mu + trasp_main + reaz_local;
        
        % Filling the sub and super diagonals
        if j > 1
            sotto(j-1) = -mu + trasp_sotto;
        end
        if j < N
            sopra(j) = -mu + trasp_sopra;
        end
    end
    
    % Construction of the tridiagonal matrix A
    A = diag(main) + diag(sopra, 1) + diag(sotto, -1);
    
    % 4. Construction of the Right-Hand Side (RHS) multiplied by h^2
    rhs = f_val * h^2;
    
    % 5. Dirichlet boundary corrections (shifting known values to the right side)
    % For node j=1, the left term (j=0) depends on alpha
    % Recalculate the specific 'sotto' coefficient that node j=1 would have had
    if eta_val(1) >= 0
        c_sotto_1 = -mu - eta_val(1) * h;
    else
        c_sotto_1 = -mu;
    end
    rhs(1) = rhs(1) - c_sotto_1 * alpha;
    
    % For node j=N, the right term (j=N+1) depends on beta
    % Recalculate the specific 'sopra' coefficient that node j=N would have had
    if eta_val(end) >= 0
        c_sopra_N = -mu;
    else
        c_sopra_N = -mu + eta_val(end) * h;
    end
    rhs(end) = rhs(end) - c_sopra_N * beta;
    
    % 6. Resolution of the linear system
    U_interni = A \ rhs;
    
    % 7. Reconstruction of the complete vector (including boundaries)
    U_com = [alpha; U_interni; beta];
end
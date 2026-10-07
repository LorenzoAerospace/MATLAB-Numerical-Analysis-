function [U_com, x] = upwind_dirichlet_1d(mu, eta, sigma, f, L, h, alpha, beta)
% UPWIND_DIRICHLET_1D  Solves a 1D steady-state advection-diffusion-reaction problem.
%
%   Solves the stationary differential equation:
%   -mu * u''(x) + eta * u'(x) + sigma * u(x) = f   in (0, L)
%   with pure Dirichlet boundary conditions: u(0) = alpha, u(L) = beta.
%   Uses a finite difference scheme with Upwind for the advection term.
%
%   INPUT:
%       mu    : Diffusion coefficient (positive constant)
%       eta   : Advection/velocity coefficient (signed constant)
%       sigma : Reaction coefficient (constant, set to 0 if absent)
%       f     : Source/forcing term (constant)
%       L     : Length of the spatial domain [0, L]
%       h     : Spatial discretization step size
%       alpha : Left boundary condition u(0) = alpha
%       beta  : Right boundary condition u(L) = beta
%
%   OUTPUT:
%       U_com : Column vector containing the complete solution (internal nodes + boundaries)
%       x     : Column vector of the spatial grid nodes

    N_intervalli = round(L/h); 
    N = N_intervalli - 1;       
    x = (0:N_intervalli)' * h;  
    
    if eta >= 0
        c_sotto = -mu - eta * h;
        c_main  =  2 * mu + eta * h + sigma * h^2; 
        c_sopra = -mu;
    else
        c_sotto = -mu;
        c_main  =  2 * mu - eta * h + sigma * h^2; 
        c_sopra = -mu + eta * h;
    end
    
    main_diag  = c_main  * ones(N, 1);
    sopra_diag = c_sopra * ones(N-1, 1);
    sotto_diag = c_sotto * ones(N-1, 1);
    A = diag(main_diag) + diag(sopra_diag, 1) + diag(sotto_diag, -1);
    
    rhs = (f * h^2) * ones(N, 1);
    
    rhs(1)   = rhs(1)   - c_sotto * alpha; 
    rhs(end) = rhs(end) - c_sopra * beta;  
    
    U_interni = A \ rhs;
    U_com = [alpha; U_interni; beta];
end
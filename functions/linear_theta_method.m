function [t, y] = linear_theta_method(A, g, t_span, y0, N, theta)
% LINEAR_THETA_METHOD  Solves the linear ODE system y' = A*y + g(t) using the theta-method.
%
%   [t, y] = LINEAR_THETA_METHOD(A, g, t_span, y0, N, theta) computes the
%   numerical solution of a system of linear ordinary differential equations 
%   using the theta-method.
%
%   INPUT PARAMETERS:
%     A      : System coefficient matrix (constant, d x d)
%     g      : Function handle for the source term g(t) -> must return a column vector
%     t_span : Vector containing initial and final times [t_start, t_end]
%     y0     : Column vector of initial conditions (d x 1)
%     N      : Number of time steps (intervals)
%     theta  : Method parameter (real constant in the interval [0, 1])
%
%   OUTPUT PARAMETERS:
%     t      : Discretized time vector ((N+1) x 1)
%     y      : Computed solution matrix ((N+1) x d)
%
%   =========================================================================
%   DOCUMENTATION OF THE THETA (\theta) PARAMETER
%   =========================================================================
%   The \theta-method is a single-step time discretization method whose 
%   nature (explicit/implicit), accuracy, and stability strictly depend on 
%   the chosen value for the parameter \theta \in [0, 1].
%
%   The three fundamental cases are:
%
%   1) \theta = 0 : FORWARD EULER METHOD (Explicit Euler)
%      ---------------------------------------------------------------------
%      * Type: Purely explicit.
%      * Accuracy order: 1st order in time O(h).
%      * Stability: Conditionally stable. The time step h must satisfy 
%        the stability condition linked to the eigenvalues of A (critical 
%        points if the system is stiff).
%      * Computational cost: Minimum per single step (no system to solve).
%
%   2) \theta = 1 : BACKWARD EULER METHOD (Implicit Euler)
%      ---------------------------------------------------------------------
%      * Type: Purely implicit.
%      * Accuracy order: 1st order in time O(h).
%      * Stability: Unconditionally stable (A-stable). Excellent for stiff 
%        or heavily damped systems. Damps high-frequency oscillations 
%        very rapidly (L-stable).
%      * Computational cost: Requires solving a linear system.
%
%   3) \theta = 0.5 : CRANK-NICOLSON METHOD (Trapezoidal Rule)
%      ---------------------------------------------------------------------
%      * Type: Implicit (Symmetric).
%      * Accuracy order: 2nd order in time O(h^2) -> Most accurate.
%      * Stability: Unconditionally stable (A-stable). However, it lacks 
%        the strong damping of Backward Euler and can produce spurious small 
%        oscillations if the step h is too large relative to the rapid 
%        dynamics of the system.
%      * Computational cost: Requires solving a linear system.
%
%   -------------------------------------------------------------------------
%   GENERAL STABILITY NOTES:
%   * For 0 <= \theta < 0.5   : The method is CONDITIONALLY stable.
%   * For 0.5 <= \theta <= 1 : The method is UNCONDITIONALLY stable.
%   =========================================================================

    % 1. Validity check for theta
    if theta < 0 || theta > 1
        error('The theta parameter must be strictly between 0 and 1.');
    end
    
    % 2. Initialization and time discretization
    t = linspace(t_span(1), t_span(2), N+1)';
    h = t(2) - t(1);
    
    % Force y0 to be a vertical column vector to avoid dimension errors
    y0 = y0(:); 
    d = length(y0);
    I = eye(d); 
    
    % Allocate solution matrix
    y = zeros(N+1, d);
    y(1, :) = y0'; 
    
    % Pre-compute the decomposition matrices for the theta-method
    M_left  = I - h * theta * A;
    M_right = I + h * (1 - theta) * A;
    
    % 3. Time integration loop
    for n = 1:N
        tn = t(n);
        tn1 = t(n+1);
        yn = y(n, :)'; % Column vector of the current state
        
        % Evaluation of g(t) at the two time nodes
        gn = g(tn);   
        gn1 = g(tn1); 
        
        % Automatic handling if the user inputs g(t) = 0 as a scalar instead of a vector
        if isscalar(gn) && gn == 0,  gn = zeros(d, 1);   end
        if isscalar(gn1) && gn1 == 0, gn1 = zeros(d, 1); end
        
        % Force g vectors to be vertical columns
        gn = gn(:); 
        gn1 = gn1(:);
        
        % Calculation of the combined known term (weighted with theta)
        rhs_vector = M_right * yn + h * theta * gn1 + h * (1 - theta) * gn;
        
        % Resolution of the linear system (efficient implicit inversion)
        y_next = M_left \ rhs_vector;
        
        % Save the result as a row
        y(n+1, :) = y_next';
    end
end
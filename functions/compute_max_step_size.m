function hmax = compute_max_step_size(A, method)
% COMPUTE_MAX_STEP_SIZE  Calculates the exact analytical h_max by solving the method's stability equation.
%
%   INPUT:
%       A      : System matrix (Jacobian)
%       method : 'EA' (Forward Euler), 'HEUN' (Heun), 'EI' (Backward Euler), 'CN' (Crank-Nicolson)
%
%   OUTPUT:
%       hmax   : The exact analytical maximum step size (scalar)

% 1. Compute the eigenvalues of the matrix
eigenvalues = eig(A);

% Implicit methods are unconditionally stable (h_max = Inf)
if strcmpi(method, 'EI') || strcmpi(method, 'CN')
    fprintf('The method %s is A-stable. h_max = Inf\n', upper(method));
    hmax = Inf;
    return;
end

% Initialize maximum step size to a very large value
hmax = Inf; 

% 2. Analytical analysis eigenvalue by eigenvalue
for i = 1:length(eigenvalues)
    lam = eigenvalues(i);
    alpha = real(lam); % Real part (\alpha)
    beta = imag(lam);  % Imaginary part (\beta)

    % Compute the squared magnitude of the eigenvalue
    mod2 = alpha^2 + beta^2; 

    switch upper(method)
        case 'EA'
            %% METHOD: FORWARD EULER
            % Exact formula from the stability circle
            h_eigenvalue = - (2 * alpha) / mod2;

        case 'HEUN'
            %% METHOD: HEUN (RK2)
            % Expansion of |1 + z + z^2/2|^2 = 1 leads to:
            c3 = (mod2^2) / 4;
            c2 = alpha * mod2;
            c1 = alpha^2 - beta^2 + mod2; 
            c0 = 2 * alpha;

            % Find the roots of the cubic polynomial: c3*h^3 + c2*h^2 + c1*h + c0 = 0
            coeff = [c3, c2, c1, c0];
            roots_poly = roots(coeff);

            % Filter for real and strictly positive roots
            real_positive_roots = roots_poly(imag(roots_poly) == 0 & roots_poly > 1e-12);

            if isempty(real_positive_roots)
                h_eigenvalue = Inf;
            else
                h_eigenvalue = min(real_positive_roots);
            end

        otherwise
            error('Unsupported method. Use ''EA'', ''HEUN'', ''EI'', or ''CN''.');
    end

    % The global maximum step size is the MINIMUM among the maximum step sizes of each eigenvalue
    if h_eigenvalue < hmax
        hmax = h_eigenvalue;
    end
end

% 3. Print final result
fprintf('=== EXACT ANALYTICAL RESULT ===\n');
fprintf('Method: %s | h_max = %.6f\n', upper(method), hmax);
fprintf('----------------------------------\n');
end
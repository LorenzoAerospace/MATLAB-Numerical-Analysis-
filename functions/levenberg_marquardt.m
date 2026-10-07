function [xvect, iter] = levenberg_marquardt(F, J, nmax, tol, x0, lambda0)
% LEVENBERG_MARQUARDT  Solves nonlinear systems/least squares using the Levenberg-Marquardt algorithm.
%
%   [XVECT, ITER] = LEVENBERG_MARQUARDT(F, J, NMAX, TOL, X0, LAMBDA0) 
%   returns the history of visited points and the number of accepted iterations.
%
%   The algorithm combines Gauss-Newton and Gradient Descent via the 
%   damping parameter LAMBDA.
%
%   INPUT:
%       F       : function handle of the system of equations
%       J       : function handle of the Jacobian matrix
%       nmax    : maximum number of iterations
%       tol     : tolerance on the step norm
%       x0      : vector of initial values
%       lambda0 : initial value of the damping parameter
%
%   OUTPUT:
%       xvect   : matrix of computed points (each row corresponds to an iteration)
%       iter    : number of accepted steps performed

% Initialization
n = numel(x0);        % Number of variables
c = 0.8;              % Lambda reduction factor
C = 2;                % Lambda increase factor
iter = 0;
x = x0;
lambda = lambda0;
xvect = x(:).';       % Store as a row for consistency in row-wise concatenation
err = tol + 1;

while (err > tol && iter < nmax)
    F_val = F(x);
    J_val = J(x);

    % Solve the normal equations (Levenberg-Marquardt step)
    % We use a weighted identity matrix to stabilize the least-squares problem
    d = (J_val' * J_val + lambda * eye(n)) \ (-J_val' * F_val);

    x_new = x + d;
    F_val_new = F(x_new);

    % Improvement test (norm comparison)
    if norm(F_val_new) < norm(F_val)
        % SUCCESS: Accept the step and reduce lambda
        x = x_new;
        lambda = c * lambda;

        % Update output parameters
        iter = iter + 1;
        err = norm(d); % Error based on step size
        xvect = [xvect; x(:).']; 
    else
        % FAILURE: Reject the step and increase lambda
        % We do not advance x, nor increment the main iteration count, but retry
        lambda = C * lambda;

        % Optional safeguard to prevent infinite loops if no improvement occurs
        if lambda > 1e10
            warning('Lambda became too large: convergence failed.');
            break;
        end
    end
end
end
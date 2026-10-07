function [t_points, y_points] = runge_kutta(f, t0, y0, t_end, h, A_butcher, b_butcher, c_butcher)
% RUNGE_KUTTA  Solves an ODE system using a generic Runge-Kutta method defined by a Butcher tableau.
%
%   [t_points, y_points] = RUNGE_KUTTA(f, t0, y0, t_end, h, A_butcher, b_butcher, c_butcher)
%
%   INPUT:
%       f         - function handle @(t, y) representing the ODE system, returning a column vector
%       t0        - initial time
%       y0        - initial state vector
%       t_end     - final integration time
%       h         - time step size
%       A_butcher - Butcher matrix (stage coefficients)
%       b_butcher - Butcher weights vector
%       c_butcher - Butcher nodes vector
%
%   OUTPUT:
%       t_points  - column vector of time steps
%       y_points  - matrix containing the computed solution at each time step (each row corresponds to a time instance)

% Ensure y0 is a column vector
y = y0(:); 
t = t0;
num_stages = length(b_butcher);
system_dim = length(y);

% Pre-allocation for efficiency (step estimation)
estimated_N = ceil((t_end - t0) ./ h) + 1;
t_points = zeros(estimated_N, 1);
y_points = zeros(estimated_N, system_dim);

% Save initial conditions
step_idx = 1;
t_points(step_idx) = t;
y_points(step_idx, :) = y';

while t < t_end
    % Handle the final step to avoid overshooting t_end
    if t + h > t_end
        h = t_end - t;
    end

    % Matrix to store k for each stage (each column is a k_i vector)
    k = zeros(system_dim, num_stages);

    % Computation of the various stages
    for i = 1:num_stages
        acc = zeros(system_dim, 1);
        for j = 1:(i-1)
            acc = acc + A_butcher(i, j) * k(:, j);
        end

        t_intermediate = t + c_butcher(i) * h;
        y_intermediate = y + h * acc;

        % Compute the derivative at this stage
        k(:, i) = f(t_intermediate, y_intermediate);
    end

    % Update of the final state for the step
    weighted_sum = zeros(system_dim, 1);
    for i = 1:num_stages
        weighted_sum = weighted_sum + b_butcher(i) * k(:, i);
    end

    y = y + h * weighted_sum;
    t = t + h;

    % Store results
    step_idx = step_idx + 1;
    t_points(step_idx) = t;
    y_points(step_idx, :) = y';
end

% Trim vectors if fewer steps were taken than estimated
t_points = t_points(1:step_idx);
y_points = y_points(1:step_idx, :);
end
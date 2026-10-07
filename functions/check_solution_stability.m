function check_solution_stability(A, h)
% CHECK_SOLUTION_STABILITY  Checks the stability of a step size h
% for a system y' = Ay + g(t) using Forward Euler, Backward Euler, Crank-Nicolson, and Heun.
%
%   INPUT:
%       A : System matrix (or scalar)
%       h : Time discretization step size (integration step)

% 1. Compute the eigenvalues of matrix A
lambda = eig(A);

% Compute z = lambda * h for each eigenvalue
Z = lambda * h;

% 2. Definition of the amplification functions R(z)
R_EA   = @(z) 1 + z;
R_EI   = @(z) 1 ./ (1 - z);
R_CN   = @(z) (1 + z/2) ./ (1 - z/2);
R_Heun = @(z) 1 + z + (z.^2)/2;
methods_name = {'Forward Euler', 'Backward Euler', 'Crank-Nicolson', 'Heun'};
functions    = {R_EA, R_EI, R_CN, R_Heun};

fprintf('===================================================\n');
fprintf(' STABILITY ANALYSIS (h = %g)\n', h);
fprintf(' Eigenvalues of A found: %s\n', mat2str(lambda, 4));
fprintf('===================================================\n\n');

% Prepare a single figure with 4 subplots (2x2)
figure('Name', ['Stability Analysis with h = ' num_to_str(h)], 'NumberTitle', 'off');

for m = 1:4
    R_fun = functions{m};
    name = methods_name{m};

    % Compute the magnitude of R(z) for all scaled eigenvalues
    R_magnitudes = abs(R_fun(Z));

    % The method is stable only if ALL eigenvalues are stable
    is_stable = all(R_magnitudes <= 1);

    % Print results to the command window
    fprintf('Method: %-18s -> ', name);
    if is_stable
        fprintf('STABLE (Max |R| = %.4f)\n', max(R_magnitudes));
    else
        fprintf('UNSTABLE ❌ (Max |R| = %.4f)\n', max(R_magnitudes));
    end

    % --- PLOTTING (SUBPLOT) ---
    subplot(2, 2, m);

    % Complex grid to map the stability region of the method
    [X, Y] = meshgrid(linspace(-3, 3, 200), linspace(-3, 3, 200));
    Z_grid = X + 1i*Y;
    R_grid = R_fun(Z_grid);

    % Stability region in light blue
    contourf(X, Y, abs(R_grid), [0 1], 'LineColor', 'none');
    colormap([0.75 0.85 0.95]); 
    hold on;

    % Cartesian axes
    xline(0, 'k--', 'LineWidth', 1);
    yline(0, 'k--', 'LineWidth', 1);

    % Plot all points z = lambda * h
    for j = 1:length(Z)
        if abs(R_fun(Z(j))) <= 1
            plot(real(Z(j)), imag(Z(j)), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
        else
            plot(real(Z(j)), imag(Z(j)), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        end
    end

    title(name);
    xlabel('Re(z)'); ylabel('Im(z)');
    grid on; axis equal;
    xlim([-3 3]); ylim([-3 3]);
end
fprintf('\nNote: Green points in the plot are STABLE, red ones are UNSTABLE.\n');
end

% Auxiliary function to convert h into a clean string for the figure title
function s = num_to_str(val)
s = num2str(val);
end
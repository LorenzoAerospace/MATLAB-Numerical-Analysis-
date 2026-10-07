function [is_admissible, R] = check_stability(method, lambda, dt)
% CHECK_STABILITY  Checks the stability of a numerical method for a given eigenvalue and time step.
%
%   [is_admissible, R] = CHECK_STABILITY(method, lambda, dt)
%
%   INPUT:
%       method : string specifying the method ('EA', 'EI', 'CN', 'Heun')
%       lambda : system eigenvalue (can be complex)
%       dt     : time discretization step size
%
%   OUTPUT:
%       is_admissible : true if stable, false otherwise
%       R             : value of the amplification factor

% Calculation of z
z = lambda * dt;

% Definition of the amplification factor R(z) based on the method
switch upper(method)
    case {'EA', 'FORWARD_EULER'}    % Forward Euler (Explicit)
        R_fun = @(z) 1 + z;
        name = 'Forward Euler';
    case {'EI', 'BACKWARD_EULER'}    % Backward Euler (Implicit)
        R_fun = @(z) 1 ./ (1 - z);
        name = 'Backward Euler';
    case {'CN', 'CRANK_NICOLSON'}    % Crank-Nicolson (Trapezoidal)
        R_fun = @(z) (1 + z/2) ./ (1 - z/2);
        name = 'Crank-Nicolson';
    case {'HEUN', 'RK2'}             % Heun (Runge-Kutta 2)
        R_fun = @(z) 1 + z + (z.^2)/2;
        name = 'Heun';
    otherwise
        error('Unrecognized method. Choose among: EA, EI, CN, Heun.');
end

% Stability calculation for the specific point
R = R_fun(z);
is_admissible = (abs(R) <= 1);

% --- PLOTTING STABILITY REGION ---
figure('Color', 'w');

% Creation of the grid in the complex plane
[X, Y] = meshgrid(linspace(-3, 3, 400), linspace(-3, 3, 400));
Z_grid = X + 1i*Y;
R_grid = R_fun(Z_grid);

% Draw the stability region (|R| <= 1)
contourf(X, Y, abs(R_grid), [0 1], 'LineColor', 'none');
colormap([0.7 0.85 1]); % Light blue color for the stable zone
hold on;

% Draw Cartesian axes
xline(0, 'k--', 'LineWidth', 1.2);
yline(0, 'k--', 'LineWidth', 1.2);

% Highlight the z point provided by the user
if is_admissible
    plot(real(z), imag(z), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    title(sprintf('%s: STABLE (\\vertR\\vert = %.3f)', name, abs(R)));
else
    plot(real(z), imag(z), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    title(sprintf('%s: UNSTABLE (\\vertR\\vert = %.3f)', name, abs(R)));
end

xlabel('Re(z)');
ylabel('Im(z)');
grid on;
axis equal;
legend('Stability Region', 'Im Axis', 'Re Axis', 'Your point z', 'Location', 'best');
hold off;
end
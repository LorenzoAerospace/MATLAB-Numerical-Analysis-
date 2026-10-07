function [fi, x_val, error_val] = lagrange_interpolation(arg1, arg2, arg3, arg4)
% LAGRANGE_INTERPOLATION  Polynomial interpolation using Lagrange's method.
%
%   [fi, x_val, error_val] = LAGRANGE_INTERPOLATION(a, b, N, f) computes the degree-N
%   interpolating polynomial for a function f over the interval [a, b].
%
%   [fi, x_val] = LAGRANGE_INTERPOLATION(x, y) computes the interpolating polynomial 
%   directly based on user-provided sample point vectors (x, y).
%
%   INPUT:
%       a, b    : Interval endpoints (if input uses a function handle).
%       N       : Degree of the polynomial (number of nodes - 1).
%       f       : Function handle (e.g., @(x) 1./(1+x.^2)).
%       x, y    : Vectors of sampled data points (if f is not available).
%
%   OUTPUT:
%       fi        : Values of the polynomial evaluated on the dense grid x_val.
%       x_val     : Fine evaluation grid (500 points) for plotting.
%       error_val : Pointwise absolute error (available only if function f is provided).

if nargin == 4
    % CASE A: Function f is provided
    a = arg1; b = arg2; N = arg3; f = arg4;
    x_nodes = linspace(a, b, N+1);
    y_nodes = f(x_nodes);
    mode = 'function';
elseif nargin == 2
    % CASE B: Only discrete data points are provided (e.g., loaded from a file)
    x_nodes = arg1;
    y_nodes = arg2;
    N = length(x_nodes) - 1;
    a = min(x_nodes); b = max(x_nodes);
    mode = 'data';
else
    error('Invalid number of input arguments.');
end

% 1. Compute the Polynomial
p = polyfit(x_nodes, y_nodes, N);

% 2. Evaluation grid for plotting
x_val = linspace(a, b, 500);
fi = polyval(p, x_val);

% 3. Error handling and plotting
figure('Color', 'w');
if strcmp(mode, 'function')
    % If f is available, compute the true error
    y_real = f(x_val);
    error_val = abs(y_real - fi);

    subplot(2, 1, 1);
    plot(x_val, y_real, 'k--', x_val, fi, 'r-', x_nodes, y_nodes, 'bo');
    title(['Function Interpolation (N=', num2str(N), ')']);
    legend('Exact', 'Polynomial', 'Nodes');
    grid on;

    subplot(2, 1, 2);
    semilogy(x_val, error_val, 'm'); 
    title('Absolute Error');
    grid on;
else
    % If only data is provided, the true error cannot be computed
    error_val = NaN; 
    plot(x_val, fi, 'r-', x_nodes, y_nodes, 'bo', 'MarkerFaceColor', 'b');
    title(['Experimental Data Interpolation (N=', num2str(N), ')']);
    legend('Lagrange Polynomial', 'Measurements (Data)');
    grid on;
end
end
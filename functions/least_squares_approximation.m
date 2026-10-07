function [fi, x_val, error_val] = least_squares_approximation(arg1, arg2, arg3, arg4, arg5)
% LEAST_SQUARES_APPROXIMATION  Polynomial least-squares approximation.
%
%   [fi, x_val, error_val] = LEAST_SQUARES_APPROXIMATION(a, b, N, m, f) approximates the 
%   function f over the interval [a, b] using N+1 nodes and a degree-m polynomial.
%
%   [fi, x_val] = LEAST_SQUARES_APPROXIMATION(x, y, m) approximates discrete data (x, y)
%   using a degree-m polynomial.
%
%   INPUT:
%       a, b    : Interval endpoints (for function handle).
%       N       : Number of intervals (N+1 points are used).
%       m       : Degree of the approximating polynomial (usually m < N).
%       f       : Handle of the function to approximate.
%       x, y    : Vectors of sampled data points.
%
%   OUTPUT:
%       fi        : Values of the polynomial evaluated on the dense grid x_val.
%       x_val     : Fine evaluation grid (500 points).
%       error_val : Pointwise absolute error (only if f is provided).

if nargin == 5
    % --- CASE A: Function Approximation ---
    a = arg1; b = arg2; N = arg3; m = arg4; f = arg5;
    x_nodes = linspace(a, b, N+1);
    y_nodes = f(x_nodes);
    mode = 'function';
elseif nargin == 3
    % --- CASE B: Experimental Data Approximation ---
    x_nodes = arg1;
    y_nodes = arg2;
    m = arg3;
    a = min(x_nodes); b = max(x_nodes);
    mode = 'data';
else
    error('Invalid inputs. Use (a, b, N, m, f) or (x, y, m).');
end

% 1. Compute the Least-Squares Polynomial
% polyfit with m < length(x)-1 automatically performs least squares
p = polyfit(x_nodes, y_nodes, m);

% 2. Generate evaluation grid
x_val = linspace(a, b, 500);
fi = polyval(p, x_val);

% 3. Plotting and Error Management
figure('Color', 'w');
if strcmp(mode, 'function')
    y_real = f(x_val);
    error_val = abs(y_real - fi);

    subplot(2, 1, 1);
    plot(x_val, y_real, 'k--', x_val, fi, 'r-', x_nodes, y_nodes, 'bo', 'LineWidth', 1.2);
    grid on; 
    legend('Exact', ['Least Squares (m=', num2str(m), ')'], 'Nodes');
    title(['Function Approximation: Nodes = ', num2str(length(x_nodes)), ', Degree m = ', num2str(m)]);

    subplot(2, 1, 2);
    semilogy(x_val, error_val, 'm'); 
    grid on;
    title('Absolute Error |f(x) - P_m(x)|');
else
    error_val = NaN; 
    plot(x_val, fi, 'r-', x_nodes, y_nodes, 'bo', 'MarkerFaceColor', 'b');
    grid on; 
    legend(['Polynomial Degree ', num2str(m)], 'Experimental Data');
    title(['Polynomial Regression (m = ', num2str(m), ')']);
end
end
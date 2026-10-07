function [dfi_coeff, dfi_val] = fast_polynomial_derivative(p, n, value)
% FAST_POLYNOMIAL_DERIVATIVE  Computes the n-th derivative of a polynomial efficiently.
%
%   [dfi_coeff, dfi_val] = FAST_POLYNOMIAL_DERIVATIVE(p, n, value)
%
%   INPUT:
%       p     - vector of polynomial coefficients (descending powers)
%       n     - order of the derivative (non-negative integer)
%       value - scalar value at which to evaluate the derivative
%
%   OUTPUT:
%       dfi_coeff - coefficients of the resulting derivative polynomial
%       dfi_val   - numerical evaluation of the derivative at the given point

% 1. ROBUSTNESS CHECKS (Input sanitization)
if isempty(p) || ~isnumeric(p)
    error('Invalid input p: must be a numeric vector of coefficients.');
end
if ~isscalar(n) || n < 0 || mod(n, 1) ~= 0
    error('Invalid input n: must be a non-negative integer.');
end
if ~isscalar(value) || ~isnumeric(value)
    error('Invalid input value: must be a numeric scalar.');
end

% Ensure p is a row vector for consistency
p_curr = p(:).'; 

% 2. LIGHTWEIGHT COMPUTATION OF THE n-TH DERIVATIVE
% If the order of derivation is greater than or equal to the number of coefficients,
% the derivative is zero.
if n >= length(p_curr)
    dfi_coeff = 0;
    dfi_val = 0;
    return;
end

% Apply polyder in a loop n times (much faster than symbolic math)
for i = 1:n
    p_curr = polyder(p_curr);
end

% 3. RESULTS
dfi_coeff = p_curr;              % Coefficients of the derivative polynomial
dfi_val = polyval(p_curr, value); % Numerical evaluation at the point
end
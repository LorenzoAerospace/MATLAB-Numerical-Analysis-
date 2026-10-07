function [succ, it] = fixed_point_iteration(x0, phi, nmax, toll)
% FIXED_POINT_ITERATION  Approximates a fixed point of the function phi using iteration.
%
%   [succ, it] = FIXED_POINT_ITERATION(x0, phi, nmax, toll)
%
%   INPUT:
%       x0   - initial iterate vector
%       phi  - function handle representing the fixed-point function phi(x)
%       nmax - maximum number of iterations
%       toll - tolerance on the stopping criterion (infinity norm of successive difference)
%
%   OUTPUT:
%       succ - matrix containing all iterate vectors column by column
%       it   - number of iterations performed

% Ensure initial format is a column vector
x0 = x0(:);
err   = 1 + toll;
it    = 0;
succ  = x0; 
xv    = x0;

while (it < nmax && err > toll)
    xn    = phi(xv);
    xn    = xn(:); % Additional column linearization check

    % Use the maximum of the absolute value (corresponds to the infinity norm)
    err   = max(abs(xn - xv)); 

    % Append the column vector xn to the succ matrix
    succ  = [succ, xn]; 
    it    = it + 1;
    xv    = xn;
end
end
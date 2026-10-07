function [xvect, it] = modified_newton(x0, nmax, toll, fun, dfun, mol)
% MODIFIED_NEWTON  Approximates a zero (root) of a function using Newton's method.
%
%   [xvect, it] = MODIFIED_NEWTON(x0, nmax, toll, fun, dfun, mol) 
%   applies the modified Newton's method for roots with a specified multiplicity (mol). 
%   The stopping criterion is based on the difference between two successive iterates.
%
%   INPUT:
%       x0        - Initial iterate
%       nmax      - Maximum number of iterations
%       toll      - Tolerance on the stopping criterion (successive iterate difference)
%       fun, dfun - Function handles containing the function and its derivative
%       mol       - Multiplicity assigned to the zero (optional, default: mol = 1)
%
%   OUTPUT:
%       xvect     - Vector containing all computed iterates (the last component is the approximated zero)
%       it        - Number of iterations performed

if (nargin == 5)
    mol = 1;
end

err = toll + 1;
it = 0;
xvect = x0;
xv = x0;
xv = xv(:);

while (it < nmax && err > toll)
    dfx = dfun(xv);
    if dfx == 0
        error('Stopped due to zero derivative (dfun = 0).');
    else
        xn = xv - mol * fun(xv) / dfx;
        xn = xn(:);
        err = max(abs(xn - xv));
        xvect = [xvect, xn];
        it = it + 1;
        xv = xn;
    end
end

if (it < nmax)
    fprintf('Converged at step k: %d\n', it);
else
    fprintf('Maximum number of steps reached: %d\n', it);
end

fprintf('Calculated root       : %-12.8f\n', xvect(end));
end
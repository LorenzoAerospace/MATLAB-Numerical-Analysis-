function [gamma_interval] = fixed_point_convergence(alpha, fun, x, beta)
% FIXED_POINT_CONVERGENCE  Determines the interval of values for the parameter beta
%                          for which the fixed-point iteration method converges.
%
%   gamma_interval = fixed_point_convergence(alpha, fun, x, beta)
%
%   INPUT:
%       alpha  - candidate fixed point (numeric value)
%       fun    - fixed-point function g(x) expressed symbolically
%       x      - symbolic variable of the function fun
%       beta   - symbolic parameter with respect to which the convergence
%                interval is to be found
%
%   OUTPUT:
%       gamma_interval - symbolic structure containing the interval of
%                        beta values satisfying the condition
%                        |g'(alpha)| < 1, along with any additional conditions
%
%   DESCRIPTION:
%       The function computes the symbolic derivative g'(x), evaluates it at
%       the fixed point alpha using subs, and solves the inequality
%           |g'(alpha)| < 1
%       with respect to the parameter beta. Symbolic variables are declared
%       as real to avoid complex or parametric solutions.

assume(x, 'real');
assume(beta, 'real');
dfun = diff(fun, x);
dfun_val = simplify(subs(dfun, x, alpha));
gamma_interval = solve(abs(dfun_val) < 1, beta, "ReturnConditions", true);
end
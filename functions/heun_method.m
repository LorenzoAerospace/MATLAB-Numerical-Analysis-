function [t, U] = heun_method(F, t0, tf, y0, h)
% HEUN_METHOD  Solves an ODE system using Heun's method (Explicit Trapezoidal / RK2).
%
%   [t, U] = HEUN_METHOD(F, t0, tf, y0, h)
%
%   INPUT:
%       F  - function handle representing the system f(t, y)
%       t0 - initial time
%       tf - final time
%       y0 - initial state vector
%       h  - time step size
%
%   OUTPUT:
%       t  - column vector of time steps
%       U  - matrix containing the computed solution at each time step (each row corresponds to a time instance)

N = round((tf - t0) / h);
t = linspace(t0, tf, N+1).';
m = length(y0);
U = zeros(N+1, m);
U(1, :) = y0(:).';

for n = 1:N
    tn = t(n);
    un = U(n, :).';

    % Heun's method (RK2) stages
    k1 = F(tn, un);
    u_tilde = un + h * k1;
    k2 = F(tn + h, u_tilde);

    % Solution update
    U(n+1, :) = (un + (h / 2) * (k1 + k2)).';
end
end
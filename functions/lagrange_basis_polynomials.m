function [phi, pol_val] = lagrange_basis_polynomials(nodes, value)
% LAGRANGE_BASIS_POLYNOMIALS  Computes the Lagrange basis polynomials.
%
%   [PHI, POL_VAL] = LAGRANGE_BASIS_POLYNOMIALS(NODES, VALUE) returns 
%   the basis (characteristic) polynomials for Lagrange interpolation 
%   built on the set of points specified in NODES.
%
%   INPUT:
%       NODES - Row or column vector containing the interpolation 
%               abscissas (x_0, x_1, ..., x_n). Nodes must be distinct.
%       VALUE - Scalar or vector of points at which to evaluate the 
%               computed polynomials.
%
%   OUTPUT:
%       PHI     - Cell array (n x 1) containing the @(x) function handles 
%                 of the basis polynomials L_k(x).
%       POL_VAL - Matrix (n x m) containing the numerical values of the 
%                 polynomials PHI evaluated at the points specified in VALUE 
%                 (where m is the length of VALUE).
%
%   DESCRIPTION:
%       The function constructs the polynomials L_k(x) such that L_k(nodes(j)) = delta_kj,
%       meaning they evaluate to 1 at the k-th node and 0 at all others.
%       It uses cell arrays and anonymous functions to ensure fast execution
%       without depending on the Symbolic Math Toolbox.
%
%   EXAMPLE:
%       nodes = [0, 1, 2];
%       val = [0.5, 1.5];
%       [p, v] = lagrange_basis_polynomials(nodes, val);
%       % To evaluate the first polynomial at x=0.5: p{1}(0.5)

% Force nodes to be a row vector for calculation consistency
nodes = nodes(:).'; 
n = length(nodes);

phi = cell(n, 1); 

for i = 1:n
    % Isolate the current node and all other nodes
    node_i = nodes(i);
    other_nodes = nodes([1:i-1, i+1:end]);

    % Define the i-th polynomial utilizing vectorized "prod".
    % arrayfun allows the function to handle vector 'x' inputs as well.
    phi{i} = @(x) arrayfun(@(xi) prod((xi - other_nodes) ./ (node_i - other_nodes)), x);
end

% Evaluation at the requested points (handles scalars, row, or column vectors)
m = length(value);
pol_val = zeros(n, m);
for k = 1:n
    pol_val(k, :) = phi{k}(value);
end
end
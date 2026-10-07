function [sse] = sum_squared_error(y_val, y_nodes)
% SUM_SQUARED_ERROR  Computes the sum of squared differences (residuals).
%
%   SSE = SUM_SQUARED_ERROR(y_val, y_nodes) returns the sum of the
%   squared differences between the vector of computed values (y_val) 
%   and the vector of reference values (y_nodes).
%
%   Inputs y_val and y_nodes must be numeric vectors of the same size. 
%   If their dimensions differ, the function returns an error.
%
%   INPUT:
%       y_val   - numeric vector of estimated/computed values
%       y_nodes - numeric vector of reference/true values
%
%   OUTPUT:
%       sse     - sum of squared errors (residuals)

% 1. Dimensional check
if numel(y_val) ~= numel(y_nodes)
    error('Inputs y_val and y_nodes must have the same number of elements.');
end

% 2. Vectorized computation
sse = sum((y_val(:) - y_nodes(:)).^2);
end
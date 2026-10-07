function [] = check_matrix(A)
% CHECK_MATRIX  Analyzes fundamental properties of a square matrix.
%
%   CHECK_MATRIX(A)
%
%   INPUT:
%       A - square matrix to analyze
%
%   The function checks the following properties:
%
%   1) Symmetry and positive definiteness (SPD):
%        - Checks if A is symmetric (A = A') 
%        - Checks if all eigenvalues are strictly positive
%
%   2) Strict row diagonal dominance:
%        |a_ii| > sum of |a_ij| for j ≠ i
%
%   3) Strict column diagonal dominance:
%        |a_jj| > sum of |a_ij| for i ≠ j
%
%   4) Non-singularity:
%        - Checks if det(A) ≠ 0
%
%   For each property, the function prints a message indicating
%   whether the matrix satisfies the condition.
n = size(A, 1);

% --- SPD Analysis ---
if isequal(A, A') && all(eig(A) > 0)
    fprintf("The matrix is SPD (Symmetric Positive Definite)\n");
else
    fprintf("The matrix is NOT SPD\n");
end

% --- Strict Row Diagonal Dominance Analysis ---
is_ddr = true;
for i = 1 : n
    % SUM OF ABSOLUTE VALUES of off-diagonal elements
    somma_ext = sum(abs(A(i, :))) - abs(A(i, i));
    if abs(A(i, i)) <= somma_ext
        is_ddr = false;
        break;
    end
end
if is_ddr
    fprintf("The matrix is strictly ROW diagonally dominant\n");
else
    fprintf("The matrix is NOT strictly ROW diagonally dominant\n");
end

% --- Strict Column Diagonal Dominance Analysis ---
is_ddc = true;
for j = 1 : n
    somma_ext = sum(abs(A(:, j))) - abs(A(j, j));
    if abs(A(j, j)) <= somma_ext
        is_ddc = false;
        break;
    end
end
if is_ddc
    fprintf("The matrix is strictly COLUMN diagonally dominant\n");
else
    fprintf("The matrix is NOT strictly COLUMN diagonally dominant\n");
end

% --- Non-singularity Analysis ---
if det(A) ~= 0
    fprintf("The matrix is NON-SINGULAR\n");
else
    fprintf("The matrix is SINGULAR\n");
end
end
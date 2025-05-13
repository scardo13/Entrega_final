function [x, Aug_final] = gaussian_elimination_simple(A, b)
% gaussian_elimination_simple Implements Simple Gaussian Elimination (no pivoting) to solve Ax=b.
% Displays augmented matrix at each stage.
% Matches the output format of the document for Method 5.
% Code elements (function name, variables, comments) are in English.
% Output text (headers, messages) is in Spanish to match the document.
%
% Inputs:
%   A - Matrix of coefficients.
%   b - Right-hand side vector.
%
% Outputs:
%   x         - Solution vector.
%   Aug_final - Final augmented matrix after elimination.

[n, m] = size(A);
if n ~= m
    error('Matrix A must be square.');
end
if size(b, 1) ~= n || size(b, 2) ~= 1
    error('Vector b must be a column vector of the same size as A.');
end

% Construct augmented matrix [A | b]
Aug = [A, b];

% Display initial augmented matrix (Etapa 0)
fprintf('Etapa 0:\n');
% Use custom print for initial augmented matrix with 6 decimal places
for row = 1:n
    for col = 1:n+1
        fprintf('%10.6f ', Aug(row, col));
    end
    fprintf('\n');
end

% Simple Gaussian elimination (no pivoting)
for k = 1:n - 1 % k is the current pivot column (from 1 to n-1)
    % Check for zero pivot Aug(k, k)
    if Aug(k, k) == 0
        error('Pivot element Aug(%d,%d) is zero. Simple Gaussian elimination fails without pivoting.', k, k);
    end

    % Perform elimination for rows below the pivot row
    for i = k+1:n % i is the current row being eliminated (from k+1 to n)
        % Calculate multiplier
        factor = Aug(i, k) / Aug(k, k);

        % Subtract factor * row k from row i
        Aug(i, k:n+1) = Aug(i, k:n+1) - factor * Aug(k, k:n+1);
    end

    % Display augmented matrix at Etapa k
    fprintf('Etapa %d:\n', k);
    % Use custom print for augmented matrix with 6 decimal places
    for row = 1:n
        for col = 1:n+1
            fprintf('%10.6f ', Aug(row, col));
        end
        fprintf('\n');
    end
end

% Backward substitution to solve for the variables
x = zeros(n, 1);
% The augmented matrix is now [U | c]
% Solve U * x = c
for i = n:-1:1 % i is the current row (from n down to 1)
    sum_ux = 0;
    for j = i+1:n % j is the column index for the sum (from i+1 to n)
        sum_ux = sum_ux + Aug(i, j) * x(j);
    end
    % Aug(i,i) is the diagonal element of the upper triangular matrix
    % Check for zero on the diagonal before division
     if Aug(i, i) == 0
         error('Error internal: Diagonal element Aug(%d,%d) is zero during backward substitution.', i, i);
     end
    x(i) = (Aug(i, n+1) - sum_ux) / Aug(i, i);
end

% Store final augmented matrix
Aug_final = Aug;

% Display final solution vector
fprintf('\nDespués de aplicar sustitución regresiva:\n');
fprintf('X:\n');
for row = 1:n
    fprintf('%10.6f\n', x(row));
end


end
function [x_final, Aug_final, variables_order_final] = gaussian_elimination_total_pivot(A, b)
% gaussian_elimination_total_pivot Implements Gaussian Elimination with Total Pivoting to solve Ax=b.
% Displays augmented matrix at each stage.
% Matches the output format of the document for Method 6.
% Code elements (function name, variables, comments) are in English.
% Output text (headers, messages) is in Spanish to match the document.
%
% Inputs:
%   A - Matrix of coefficients.
%   b - Right-hand side vector.
%
% Outputs:
%   x_final             - Solution vector.
%   Aug_final           - Final augmented matrix after elimination.
%   variables_order_final - Final order of variables after column swaps.

[n, m] = size(A);
if n ~= m
    error('Matrix A must be square.');
end
if size(b, 1) ~= n || size(b, 2) ~= 1
    error('Vector b must be a column vector of the same size as A.');
end

% Construct augmented matrix [A | b]
Aug = [A, b];

% Initialize vector to track original variable order (1 to n)
variables_order = 1:n;

% Display initial augmented matrix (Etapa 0)
fprintf('Etapa 0:\n');
% Use custom print for initial augmented matrix with 6 decimal places
for row = 1:n
    for col = 1:n+1
        fprintf('%10.6f ', Aug(row, col));
    end
    fprintf('\n');
end

% Gaussian elimination with total pivoting
for k = 1:n % k is the current pivot column/row (from 1 to n)
    % Find the pivot element (largest absolute value in the submatrix Aug(k:n, k:n))
    submatrix = Aug(k:n, k:n);
    [max_val_in_sub, linear_idx_in_sub] = max(abs(submatrix(:))); % Find max value and its linear index in the submatrix
    [row_idx_sub, col_idx_sub] = ind2sub(size(submatrix), linear_idx_in_sub); % Convert linear index to row/col in submatrix

    % Convert submatrix indices to global indices in the augmented matrix
    pivot_row = row_idx_sub + k - 1;
    pivot_col = col_idx_sub + k - 1;

    % Check if the pivot element is zero (matrix is singular or close to singular)
    if max_val_in_sub == 0
        error('Matrix is singular or close to singular. Gaussian elimination with total pivoting fails at stage %d.', k);
    end

    % Swap rows k and pivot_row (if necessary)
    if pivot_row ~= k
        Aug([k, pivot_row], :) = Aug([pivot_row, k], :);
    end

    % Swap columns k and pivot_col (if necessary)
    if pivot_col ~= k
        Aug(:, [k, pivot_col]) = Aug(:, [pivot_col, k]);
        % Swap the corresponding entries in the variables order vector
        variables_order([k, pivot_col]) = variables_order([pivot_col, k]);
    end

    % Perform elimination for rows below the pivot row
    for i = k+1:n
        % Calculate multiplier
        % Check if the pivot element (Aug(k,k)) is zero before division
         if Aug(k, k) == 0 % This should be caught by max_val_in_sub check, but as a safeguard
             error('Internal error: Pivot element is zero during elimination at stage %d.', k);
         end
        factor = Aug(i, k) / Aug(k, k); % Aug(k,k) is the pivot element after swaps

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

% Backward substitution to solve for the permuted variables
x_permuted = zeros(n, 1);
% The augmented matrix is now in row echelon form [U_permuted | c_permuted], but with columns permuted
% Solve U_permuted * x_permuted = c_permuted
for i = n:-1:1
    sum_ux = 0;
    for j = i+1:n
        sum_ux = sum_ux + Aug(i, j) * x_permuted(j);
    end
    % Aug(i,i) is the diagonal element of the upper triangular matrix (after column swaps)
    % Check for zero on the diagonal before division
     if Aug(i, i) == 0
         error('Error internal: Diagonal element Aug(%d,%d) is zero during backward substitution.', i, i);
     end
    x_permuted(i) = (Aug(i, n+1) - sum_ux) / Aug(i, i);
end

% Rearrange the solution vector back to the original variable order
x_final = zeros(n, 1);
for i = 1:n
    % The variable that was originally at index variables_order(i) is now at index i in x_permuted
    % So, x_final at the original index should get the value from x_permuted at the current index i
    x_final(variables_order(i)) = x_permuted(i);
end

% Store final outputs
Aug_final = Aug;
variables_order_final = variables_order;

% Display final solution vector
fprintf('\nDespués de aplicar sustitución regresiva:\n');
fprintf('X:\n');
for row = 1:n
    fprintf('%10.6f\n', x_final(row));
end


end
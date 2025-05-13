function [root, iter_out, error_hist] = secant_method(func, x0, x1, tol, nmax)
% secant_method Implements the Secant method for finding a root.
% Uses absolute error |xn - xn-1| as stopping criterion.
% Displays a table of iterations matching the document format.
% Code elements (function name, variables, comments) are in English.
% Output text (headers, messages) is in Spanish to match the document.
%
% Inputs:
%   func  - Function handle of the function f(x).
%   x0, x1 - Two initial guesses for the root.
%   tol   - Tolerance for the stopping criterion (absolute error |xn - xn-1|).
%   nmax  - Maximum number of iterations.
%
% Outputs:
%   root       - Approximation of the root found.
%   iter_out   - Number of iterations performed (last iteration number displayed).
%   error_hist - History of errors in each iteration.

% Prepare the table header (matching the document format)
fprintf('iter |         xi        |  f(xi)  |    E    |\n');
fprintf('-----|-------------------|---------|---------|\n');

% Initialize variables with the two initial guesses
xn_minus_1 = x0; % x_{n-1} (initial guess 1)
xn = x1;         % x_n (initial guess 2)
error = Inf;     % Initialize error to a value greater than tol
error_hist = []; % Initialize error history

% Print rows for initial guesses (iter 0 and 1) as per document format
fprintf('%4d | %17.10f | %+8.1e |         |\n', 0, xn_minus_1, func(xn_minus_1));
fprintf('%4d | %17.10f | %+8.1e |         |\n', 1, xn, func(xn));

% Iteration loop (calculate x_2, x_3, ...)
% The loop will perform iterations from 2 up to nmax.
% Use iter as the loop variable and the table iteration number.

for iter = 2:nmax % iter is the current iteration number displayed in the table (starts from 2)

    % Evaluate function at current and previous approximations
    f_xn = func(xn);
    f_xn_minus_1 = func(xn_minus_1);

    % Check for zero denominator (f(xn) - f(xn-1)) in the Secant formula
    denominator = f_xn - f_xn_minus_1;
    if abs(denominator) < eps
        warning('El denominador en la fórmula de la Secante es cercano a cero en la iteración %d. El método podría tener problemas.', iter);
         if denominator == 0
             error('El denominador en la fórmula de la Secante es cero en la iteración %d.', iter);
         end
    end

    % Calculate the next approximation (x_next) using Secant formula
    x_next = xn - f_xn * (xn - xn_minus_1) / denominator;

    % Calculate absolute error |xn - xn-1| (where xn is x_next, xn-1 is xn from current iter)
    error = abs(x_next - xn);
    error_hist = [error_hist, error]; %#ok<AGROW> % Store error history

    % Evaluate function at the next approximation (for printing)
    f_x_next = func(x_next);

    % Print table row
    fprintf('%4d | %17.10f | %+8.1e | %+8.1e |\n', iter, x_next, f_x_next, error);

    % Update approximations for the next iteration
    xn_minus_1 = xn;
    xn = x_next;

    % Check stopping criterion based on absolute error
    if error < tol
        root = xn; % The current approximation (x_next, now stored in xn) is the root
        % Removed convergence message fprintf here
        iter_out = iter; % Set output iteration count
        break; % Exit the loop due to convergence
    end

    % Optional: Check if f(xn) (which is f_x_next) is effectively zero (found the root exactly or very close)
     if abs(f_x_next) < eps
         root = xn;
         % Removed convergence message fprintf here
         iter_out = iter; % Set output iteration count
         break; % Exit the loop
     end

end % End of for loop

% After the loop finishes (either by break or reaching nmax)

% Determine the final iteration count for the output variable iter_out
% If the loop broke, iter_out was set inside the if block.
% If the loop finished naturally, iter reached nmax.
if iter == nmax && error >= tol % Check if loop finished naturally without convergence (iter reaches nmax)
     iter_out = nmax; % Set output iteration count to nmax
     warning('El método no convergió después de %d iteraciones con la tolerancia de error absoluto especificada.', nmax);
end
% If loop broke, iter_out was already set inside the if block.

root = xn; % The last calculated approximation is the root

% Print the final root approximation matching the document's format
% This line will now be reached after a break or after the loop finishes naturally.
fprintf('\nSe encontró una aproximación de la raiz en %.15f\n', xn);

end
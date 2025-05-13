function [root, iter, error_hist] = newton_method(func, dfunc, x0, tol, nmax)
% newton_method Implements the Newton-Raphson method for finding a root.
% Uses absolute error |xn - xn-1| as stopping criterion.
% Displays a table of iterations matching the document format.
% Code elements (function name, variables, comments) are in English.
% Output text (headers, messages) is in Spanish to match the document.
%
% Inputs:
%   func  - Function handle of the function f(x).
%   dfunc - Function handle of the derivative of f(x), f'(x).
%   x0    - Initial guess for the root.
%   tol   - Tolerance for the stopping criterion (absolute error |xn - xn-1|).
%   nmax  - Maximum number of iterations.
%
% Outputs:
%   root       - Approximation of the root found.
%   iter       - Number of iterations performed.
%   error_hist - History of errors in each iteration.

% Prepare the table header (matching the document format)
fprintf('iter |         xi        |  f(xi)  |    E    |\n');
fprintf('-----|-------------------|---------|---------|\n');

% Initialize variables
iter = 0;
x_current = x0;
error = Inf; % Initialize error to a value greater than tol to enter the loop
error_hist = []; % Initialize error history

% Print row for iteration 0 (initial guess)
f_x0 = func(x0);
fprintf('%4d | %17.10f | %+8.1e |         |\n', iter, x0, f_x0); % No error for iter 0

% Iteration loop
while iter < nmax
    iter = iter + 1; % Increment iteration counter

    % Evaluate function and its derivative at the current approximation
    f_x = func(x_current);
    df_x = dfunc(x_current);

    % Check for zero derivative (potential division by zero)
    if abs(df_x) < eps % Using machine epsilon for near-zero check
        warning('La derivada es cero o muy cercana a cero en la iteración %d. El método de Newton-Raphson falla.', iter);
        root = x_current; % The last approximation is the best found before failure
        return; % Exit the function due to failure
    end

    % Calculate the next approximation using Newton-Raphson formula: x_next = x_current - f(x_current) / f'(x_current)
    x_next = x_current - f_x / df_x;

    % Calculate absolute error |xn - xn-1| (where xn is x_next, xn-1 is x_current)
    error = abs(x_next - x_current);
    error_hist = [error_hist, error]; %#ok<AGROW> % Store error history

    % Evaluate function at the next approximation (for printing in the table)
    f_x_next = func(x_next);

    % Print table row
    fprintf('%4d | %17.10f | %+8.1e | %+8.1e |\n', iter, x_next, f_x_next, error);

    % Update x_current for the next iteration
    x_current = x_next;

    % Check stopping criterion based on absolute error
    if error < tol
        root = x_current; % The current approximation is the root
        % Removed convergence message fprintf here
        break; % Exit the loop due to convergence
    end

    % Optional: Check if f(x_current) is effectively zero (found the root exactly or very close)
     if abs(f_x_next) < eps
         root = x_current;
         % Removed convergence message fprintf here
         break; % Exit the loop
     end

end % End of while loop

% If the loop finishes because the maximum number of iterations was reached AND error tolerance was not met
if iter == nmax && error >= tol
    root = x_current; % The last approximation is the best found
    warning('El método no convergió después de %d iteraciones con la tolerancia de error absoluto especificada.', nmax);
end
% Note: If the method converged or failed due to zero derivative, the loop was exited early.

% Print the final root approximation matching the document's format
fprintf('\nSe encontró una aproximación de la raiz en %.15f\n', root);

end
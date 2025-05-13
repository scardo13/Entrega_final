function [root, iter, error_hist] = bisection_method(func, a, b, tol, nmax)
% bisection_method Implements the Bisection method for finding a root of func in [a, b].
% Uses absolute error |xn - xn-1| as stopping criterion.
% Displays a table of iterations matching the document format.
%
% Inputs:
%   func  - Function handle of the function f(x).
%   a, b  - Endpoints of the initial interval [a, b].
%   tol   - Tolerance for the stopping criterion (absolute error |xn - xn-1|).
%   nmax  - Maximum number of iterations.
%
% Outputs:
%   root       - Approximation of the root found.
%   iter       - Number of iterations performed.
%   error_hist - History of errors in each iteration.

% Input validation: Check if the function has opposite signs at the endpoints
if func(a) * func(b) >= 0
    error('La función debe tener signos opuestos en los extremos a y b para asegurar la existencia de una raíz.');
end

% Prepare the table header (matching the document format)
fprintf('iter |         a         |         xm        |         b         |  f(xm)  |    E    |\n');
fprintf('-----|-------------------|-------------------|-------------------|---------|---------|\n');

% Initialize variables
iter = 0;
error = Inf; % Initialize error to a value greater than tol to enter the loop
xm_prev = a; % Initialize previous midpoint with 'a' for the first error calculation |x1 - x0| where x0 = a
root = a; % Initialize root, will be updated in the loop
error_hist = []; % Initialize error history vector

% Iteration loop: Continue until max iterations reached or error tolerance met
while iter < nmax && error >= tol
    iter = iter + 1; % Increment iteration counter

    % Calculate midpoint
    xm = a + (b - a) / 2; % Calculate the midpoint of the current interval [a, b]

    % Evaluate function at midpoint
    f_xm = func(xm);

    % Calculate absolute error |xn - xn-1|
    % For iter 1, xn is xm and xn-1 is xm_prev (which is 'a')
    % For iter > 1, xn is the current xm and xn-1 is the xm from the previous iteration
    error = abs(xm - xm_prev);
    error_hist = [error_hist, error]; %#ok<AGROW> % Store error history

    % Update previous midpoint for the next iteration
    xm_prev = xm;

    % Print table row (matching the document's formatting)
    % iter: integer
    % a, xm, b: fixed point with 10 decimal places
    % f(xm), E: scientific notation with 1 digit after decimal point (e.g., -2.9e-01)
    fprintf('%4d | %17.10f | %17.10f | %17.10f | %+8.1e | %+8.1e |\n', iter, a, xm, b, f_xm, error);

    % Check stopping criterion based on absolute error
    if error < tol
        root = xm; % The current midpoint is the root approximation
        fprintf('\nConvergencia alcanzada por error absoluto en %d iteraciones.\n', iter);
        % Optional: Also check if f(xm) is very close to zero, as an alternative convergence sign
        if abs(f_xm) < eps % eps is machine epsilon
             fprintf(' (y el valor de la función está cercano a cero)');
        end
        fprintf('\n');
        % The final result message is printed after the loop, so no need for a return here
    end

    % Update the interval for the next iteration
    if func(a) * f_xm < 0
        % The root is in the left subinterval [a, xm]
        b = xm;
    else
        % The root is in the right subinterval [xm, b] (or xm is the root if f_xm = 0)
        a = xm;
    end

    % Optional: Check if f(xm) is effectively zero (found the root exactly or very close)
    if abs(f_xm) < eps % Using machine epsilon for near-zero check
         root = xm;
         fprintf('\nConvergencia alcanzada porque el valor de la función está cercano a cero en %d iteraciones.\n', iter);
         % The final result message is printed after the loop
         break; % Exit the loop if f(xm) is close to zero
    end
end % End of while loop

% If the loop finishes because the maximum number of iterations was reached AND error tolerance was not met
if iter == nmax && error >= tol
    root = xm; % The last calculated midpoint is the best approximation
    warning('El método no convergió después de %d iteraciones con la tolerancia de error absoluto especificada.', nmax);
end
% Note: If convergence was met by error or f(xm) close to zero, the loop was exited early,
% and 'root' was set to the converged value.

% Print the final root approximation matching the document's format
% This line will execute whether convergence was met or max iterations were reached.
fprintf('\nSe encontró una aproximación de la raiz en %.15f\n', root);

end % End of function
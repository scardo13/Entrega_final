function [root, iter, error_hist] = false_position_method(func, a, b, tol, nmax)
% false_position_method Implements the False Position method for finding a root of func in [a, b].
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
f_a_init = func(a);
f_b_init = func(b);
if f_a_init * f_b_init >= 0
    error('La función debe tener signos opuestos en los extremos a y b para asegurar la existencia de una raíz.');
end

% Prepare the table header (matching the document format)
fprintf('iter |         a         |         xm        |         b         |  f(xm)  |    E    |\n');
fprintf('-----|-------------------|-------------------|-------------------|---------|---------|\n');

% Initialize variables
iter = 0;
error = Inf; % Initialize error to a value greater than tol to enter the loop
xm_prev = a; % Initialize previous xm with 'a' for the first error calculation |x1 - x0| where x0 = a
root = a; % Initialize root
error_hist = []; % Initialize error history

% Iteration loop
while iter < nmax && error >= tol
    iter = iter + 1;

    % Evaluate function at current endpoints
    f_a = func(a);
    f_b = func(b);

    % Calculate the next approximation (xm) using False Position formula
    % Formula: xm = b - f(b) * (b - a) / (f(b) - f(a))
    % Check for division by zero f(b) - f(a)
    denominator = f_b - f_a;
    if abs(denominator) < eps % Check if denominator is close to zero
        warning('El denominador en la fórmula de Regla Falsa es cercano a cero. El método podría tener problemas de convergencia o el intervalo no es adecuado.');
        % If the denominator is exactly zero, the method fails
        if denominator == 0
             error('El denominador en la fórmula de Regla Falsa es cero.');
        end
    end
    xm = b - f_b * (b - a) / denominator;

    % Evaluate function at xm
    f_xm = func(xm);

    % Calculate absolute error |xn - xn-1|
     if iter > 1
        error = abs(xm - xm_prev);
    else % For the first iteration, error is |x1 - x0|, where x1=xm, x0=a
        error = abs(xm - a); % Using 'a' as the previous point for iter 1 error
    end
    error_hist = [error_hist, error]; %#ok<AGROW>

    % Update previous xm for the next iteration
    xm_prev = xm;

    % Print table row
    fprintf('%4d | %17.10f | %17.10f | %17.10f | %+8.1e | %+8.1e |\n', iter, a, xm, b, f_xm, error);

    % Check stopping criterion based on absolute error
    if error < tol
        root = xm;
        fprintf('\nConvergencia alcanzada por error absoluto en %d iteraciones.\n', iter);
         if abs(f_xm) < eps
             fprintf(' (y el valor de la función está cercano a cero)');
        end
        fprintf('\n');
        break; % Exit loop
    end

    % Update interval: Check where the root lies based on signs
    if f_a * f_xm < 0
        b = xm; % Root is in [a, xm]
    else % f_xm * f_b < 0 or f_xm = 0
        a = xm; % Root is in [xm, b]
    end

    % Optional: Check if f(xm) is effectively zero
    if abs(f_xm) < eps
         root = xm;
         fprintf('\nConvergencia alcanzada porque el valor de la función está cercano a cero en %d iteraciones.\n', iter);
         break; % Exit loop
    end
end % End of while loop

% If loop finishes because max iterations reached AND error tolerance was not met
if iter == nmax && error >= tol
    root = xm; % Last calculated xm
    warning('El método no convergió después de %d iteraciones con la tolerancia de error absoluto especificada.', nmax);
end

% Print final root approximation
fprintf('\nSe encontró una aproximación de la raiz en %.15f\n', root);

end
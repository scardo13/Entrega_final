function [Li_polynomials_coefs, final_polynomial_expr_str] = lagrange_interpolation(x_values, y_values)
% lagrange_interpolation Implements Lagrange interpolation.
% Calculates and displays individual Lagrange basis polynomials L_i(x)
% and the final interpolating polynomial P(x) = sum(y_i * L_i(x)).
% Matches the output format of the document for Method 11.
%
% Entradas:
%   x_values - Vector de valores de x.
%   y_values - Vector de valores de y correspondientes.
%
% Salidas:
%   Li_polynomials_coefs  - Cell array of polynomial coefficient vectors for each L_i(x).
%   final_polynomial_expr_str - String representing the final interpolating polynomial expression.

n = length(x_values); % Number of points

if length(y_values) ~= n
    error('Los vectores x y y deben tener la misma longitud.');
end

% Inicializar variables de salida
Li_polynomials_coefs = cell(n, 1);
final_polynomial_expr_str = '';

fprintf('Resultados:\n\n'); % Encabezado "Resultados:"

% Calculate and display individual Lagrange basis polynomials L_i(x)
fprintf('Polinomios interpolantes de Lagrange:\n'); % Encabezado

Li_strings = cell(n, 1); % To store polynomial strings for L_i

for i = 1:n % Loop through each data point (from 1 to n)
    % Roots of the numerator polynomial for L_i(x) are all x_j where j != i
    % MATLAB indexing is 1-based, so x_values(i) is the current point x_i
    % The roots are x_values(1) ... x_values(i-1) and x_values(i+1) ... x_values(n)
    if i == 1
        roots_num = x_values(2:n);
    elseif i == n
        roots_num = x_values(1:n-1);
    else
        roots_num = [x_values(1:i-1), x_values(i+1:n)];
    end


    % Coefficients of the numerator polynomial (using poly(roots))
    % poly(r) returns coefficients of the polynomial with roots r
    coef_num = poly(roots_num);

    % Denominator for L_i(x) = prod(x_i - x_j) for j != i
    denom_values = x_values(i) - roots_num;
    denom = prod(denom_values);

    % Handle potential division by zero (duplicate x_values)
    if denom == 0
        error('Puntos x duplicados detectados. No se pueden calcular los polinomios de Lagrange.');
    end

    % Coefficients of L_i(x)
    coef_Li = coef_num / denom;
    Li_polynomials_coefs{i} = coef_Li; % Store the coefficients

    % Convert L_i(x) coefficients to a formatted polynomial string
    Li_str = poly_coef_to_string(coef_Li, n-1); % Use helper function
    Li_strings{i} = Li_str; % Store the string

    % Display L_i(x) polynomial string
    fprintf('%s //L%d\n', Li_str, i-1); % Document uses L0, L1, etc. (0-based indexing)
end

% Construct the final interpolating polynomial expression string P(x) = sum(y_i * L_i(x))
fprintf('\nPolinomio:\n'); % Encabezado

final_poly_str = '';
first_term_added = false; % Flag to handle the first term's sign

for i = 1:n % Loop through each data point (from 1 to n)
    yi = y_values(i); % y_i value
    Li_name = ['L' num2str(i-1)]; % Name like L0, L1, L2, L3

    % Skip term if y_i is zero (or very close to zero)
    if abs(yi) < 1e-10
        continue;
    end

    % Determine the term string (y_i * L_i_name)
    term_str = '';
    if abs(abs(yi) - 1) < 1e-10 % Handle y_i = 1 or y_i = -1
        if yi > 0
            term_str = Li_name; % e.g., L3
        else % yi < 0
            term_str = ['-' Li_name]; % e.g., -L0 if y0=-1
        end
    else % Handle y_i is not 1 or -1
         term_str = [sprintf('%.6f', yi) '*' Li_name]; % e.g., 15.500000*L0
    end

    % Add sign and term to the final polynomial string
    if ~first_term_added
        % This is the first non-zero term
        final_poly_str = term_str; % Add the term string directly (it includes its own sign if negative)
        first_term_added = true;
    else
        % Subsequent terms
        if yi > 0
            final_poly_str = [final_poly_str '+' term_str]; % Add '+' before the term
        else
            % The '-' is already included in term_str if yi < 0, just concatenate
            final_poly_str = [final_poly_str term_str];
        end
    end
end

final_polynomial_expr_str = final_poly_str; % Store the final string

fprintf('%s\n', final_polynomial_expr_str);


end

% Helper function to convert polynomial coefficients to a formatted string
% coefs: vector of coefficients [c_k, c_{k-1}, ..., c_1, c_0] for degree k
% degree: the degree of the polynomial
function poly_str = poly_coef_to_string(coefs, degree)
    poly_str = '';
    n_coefs = length(coefs); % Number of coefficients (degree + 1)

    for j = 1:n_coefs
        coef = coefs(j);
        power = degree - (j - 1); % Current power of x

        % Ignore terms with coefficients very small (close to zero), except for the constant term
        if abs(coef) < 1e-10 && power ~= 0
            continue;
        end

        % Add sign (if not the very first term and coefficient is positive)
        if coef > 0 && ~isempty(poly_str) && ~endsWith(poly_str, '-') % Add '+' only if it's not the first term and not immediately after a '-'
             poly_str = [poly_str '+'];
        end
        if coef < 0
            poly_str = [poly_str '-'];
            coef = abs(coef); % Use absolute value after adding sign
        end

        % Add coefficient numerical value (only if not 1 or -1, except for constant term 0)
        % Document shows coefficient even if 1, let's follow that, but handle the constant 0
        if abs(coef) > 1e-10 || power == 0 % Always print constant term if it's non-zero, or if it's the zero constant term
             poly_str = [poly_str sprintf('%.6f', coef)];
        elseif abs(coef) <= 1e-10 && power == 0 % Explicitly print 0.000000 for zero constant term if it's the only term
             if isempty(poly_str) || endsWith(poly_str, '+') || endsWith(poly_str, '-') % Only add 0.000000 if needed to form a term
                  poly_str = [poly_str sprintf('%.6f', 0.0)];
             end
        end


        % Add x term with power
        if power > 1
            poly_str = [poly_str 'x^' num2str(power)];
        elseif power == 1
            poly_str = [poly_str 'x'];
        end
        % If power is 0 (constant term), the coefficient is already added.
    end

    % Handle case where all coefficients were effectively zero
     if isempty(poly_str)
         poly_str = sprintf('%.6f', 0.0);
     end

    % Final cleanup: ensure correct spacing and signs
    % (This can be complex to match exactly, the current logic aims for mathematical correctness close to the format)


end
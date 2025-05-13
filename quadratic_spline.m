function [coefficients, spline_strings] = quadratic_spline(x_values, y_values)
% quadratic_spline Implements quadratic spline interpolation with S_1''(x_1)=0.
% Calculates and displays coefficients and polynomial strings for each quadratic segment.
% Matches the output format of the document for Method 13.
%
% Entradas:
%   x_values - Vector de valores de x.
%   y_values - Vector de valores de y correspondientes.
%
% Salidas:
%   coefficients   - Matrix of coefficients [a_k, b_k, c_k] for each segment.
%   spline_strings - Cell array of polynomial strings for each segment.

n = length(x_values); % Number of points
if length(y_values) ~= n
    error('Los vectores x y y deben tener la misma longitud.');
end

num_segments = n - 1; % Number of quadratic segments

if num_segments <= 0
    error('Se necesitan al menos 2 puntos para crear un trazador cuadrático.');
end

% Calculate slopes m_k and segment lengths h_k
m_values = zeros(num_segments, 1);
h_values = zeros(num_segments, 1);
for k = 1:num_segments
    if (x_values(k+1) - x_values(k)) == 0
         error('Puntos x duplicados detectados entre x(%d) y x(%d). No se puede calcular el trazador cuadrático.', k, k+1);
    end
    m_values(k) = (y_values(k+1) - y_values(k)) / (x_values(k+1) - x_values(k));
    h_values(k) = x_values(k+1) - x_values(k);
end

% Initialize coefficient arrays
a_coefs = zeros(num_segments, 1);
b_coefs = zeros(num_segments, 1);
c_coefs = zeros(num_segments, 1);

% --- Calculate 'a' coefficients ---
% Extra condition: a_1 = 0
a_coefs(1) = 0;

% Recurrence relation for a_{k+1}: a_{k+1} = (m_{k+1} - m_k - a_k h_k) / h_{k+1} for k=1, ..., n-2
for k = 1:num_segments - 1
    a_coefs(k+1) = (m_values(k+1) - m_values(k) - a_coefs(k) * h_values(k)) / h_values(k+1);
end

% --- Calculate 'b' coefficients ---
% b_k = m_k - a_k * (x_{k+1} + x_k)
for k = 1:num_segments
    b_coefs(k) = m_values(k) - a_coefs(k) * (x_values(k+1) + x_values(k));
end

% --- Calculate 'c' coefficients ---
% c_k = y_k - a_k x_k^2 - b_k x_k
for k = 1:num_segments
    c_coefs(k) = y_values(k) - a_coefs(k) * x_values(k)^2 - b_coefs(k) * x_values(k);
end

% Combine coefficients into a matrix
coefficients = [a_coefs, b_coefs, c_coefs];

% Initialize cell array for spline strings
spline_strings = cell(num_segments, 1);

% Construct polynomial strings for each segment
for k = 1:num_segments
    spline_strings{k} = quadratic_poly_to_string(a_coefs(k), b_coefs(k), c_coefs(k)); % Use helper function
end

% --- Display results in the order of the document ---

fprintf('Resultados:\n\n'); % Encabezado "Resultados:"

fprintf('Coeficientes de los trazadores:\n'); % Encabezado
% Print the coefficients matrix row by row
for row = 1:num_segments
    fprintf('%10.6f %10.6f %10.6f\n', coefficients(row, 1), coefficients(row, 2), coefficients(row, 3)); % Print a_k, b_k, c_k
end

fprintf('\nTrazadores:\n'); % Encabezado
% Print each spline polynomial string
for k = 1:num_segments
    fprintf('%s\n', spline_strings{k});
end


end

% Helper function to convert polynomial coefficients [a, b, c] to a formatted string ax^2 + bx + c
function poly_str = quadratic_poly_to_string(a, b, c)
    poly_str = '';

    % Add 'ax^2' term (always include as per document format)
    poly_str = [poly_str sprintf('%.6f', a) 'x^2']; %#ok<AGROW>

    % Add 'bx' term
    if abs(b) >= 1e-10 % If b is not effectively zero
        if b > 0 % Add '+' sign if b is positive and it's not the first term (ax^2 is always first)
             poly_str = [poly_str '+']; %#ok<AGROW>
        end
        poly_str = [poly_str sprintf('%.6f', b) 'x']; %#ok<AGROW> % sprintf includes '-' sign if b is negative
    end

    % Add 'c' term (constant)
    if abs(c) >= 1e-10 % If c is not effectively zero
        % Add '+' sign if c is positive and there were previous terms (ax^2 or bx)
        if c > 0 && (abs(a) >= 1e-10 || abs(b) >= 1e-10)
             poly_str = [poly_str '+']; %#ok<AGROW>
        end
         poly_str = [poly_str sprintf('%.6f', c)]; %#ok<AGROW> % sprintf includes '-' sign if c is negative
    elseif abs(a) < 1e-10 && abs(b) < 1e-10 % If all a and b are zero, the constant term is the only one
         poly_str = sprintf('%.6f', 0.0); % In this case, print the zero constant term
    end
     % Handle case where only constant term exists and it's positive (no leading +)
     if abs(a) < 1e-10 && abs(b) < 1e-10 && c > 0 && startsWith(poly_str, '+')
          poly_str = poly_str(2:end);
     end
    % Handle case where only bx term exists and b > 0 (no leading +)
     if abs(a) < 1e-10 && abs(b) >= 1e-10 && b > 0 && startsWith(poly_str, '+')
          poly_str = poly_str(2:end);
     end


end
function [coefficients, spline_strings] = cubic_spline(x_values, y_values)
% cubic_spline Implements natural cubic spline interpolation.
% Calculates and displays coefficients and polynomial strings for each cubic segment.
% Results may not match document's coefficients due to likely use of a different spline type or algorithm in the document.
%
% Entradas:
%   x_values - Vector de valores de x.
%   y_values - Vector de valores de y correspondientes.
%
% Salidas:
%   coefficients   - Matrix of coefficients [A_k, B_k, C_k, D_k] for each segment (in x^3, x^2, x, 1 basis).
%   spline_strings - Cell array of polynomial strings for each segment.

n = length(x_values); % Number of points
if length(y_values) ~= n
    error('Los vectores x y y deben tener la misma longitud.');
end

num_segments = n - 1; % Number of cubic segments

if num_segments <= 1 % Need at least 3 points for a cubic spline (2 segments)
    error('Se necesitan al menos 3 puntos para crear un trazador cúbico (n >= 3).');
end

% Calculate segment lengths h_k and slopes m_k
h_values = zeros(num_segments, 1);
m_values = zeros(num_segments, 1);
for k = 1:num_segments
    if (x_values(k+1) - x_values(k)) == 0
         error('Puntos x duplicados detectados entre x(%d) y x(%d). No se puede calcular el trazador cúbico.', k, k+1);
    end
    h_values(k) = x_values(k+1) - x_values(k);
    m_values(k) = (y_values(k+1) - y_values(k)) / h_values(k);
end

% --- Solve for second derivatives at knots (z_k) ---
% Natural spline boundary conditions: z_1 = 0, z_n = 0
% We solve for z_k for k = 1, ..., n

z_values = zeros(n, 1); % z_k are the second derivatives at the knots x_k

% The system for z_k for k=2, ..., n-1 is tridiagonal:
% h_{k-1} z_{k-1} + 2(h_{k-1} + h_k) z_k + h_k z_{k+1} = 6 * ((y_k - y_{k-1}) / h_{k-1} - (y_{k+1} - y_k) / h_k)
% Let's use m_values = (y_{k+1}-y_k)/h_k
% h_{k-1} z_{k-1} + 2(h_{k-1} + h_k) z_k + h_k z_{k+1} = 6 * (m_{k-1} - m_k)  for k=2, ..., n-1

% Construct the tridiagonal system matrix and right-hand side vector
% The system is for z_k for k = 2, ..., n-1 (system size n-2)
system_size = n - 2;
A_sys = zeros(system_size);
b_sys = zeros(system_size, 1);

% Indices in the system matrix correspond to knots k = 2, ..., n-1
% Row i of A_sys and b_sys corresponds to knot k = i + 1

for i = 1:system_size % i refers to the row/equation in the system (1 to n-2)
    k = i + 1; % k refers to the knot index (2 to n-1)

    % Diagonal element (coefficient of z_k): 2 * (h_{k-1} + h_k)
    A_sys(i, i) = 2 * (h_values(k-1) + h_values(k));

    % Subdiagonal element (coefficient of z_{k-1}): h_{k-1}
    if i > 1
        A_sys(i, i - 1) = h_values(k-1);
    end

    % Superdiagonal element (coefficient of z_{k+1}): h_k
    if i < system_size
        A_sys(i, i + 1) = h_values(k);
    end

    % Right side (b_sys(i)): 6 * (m_k - m_{k-1})
    b_sys(i) = 6 * (m_values(k) - m_values(k-1));
end

% Solve the tridiagonal system for z_2, ..., z_{n-1}
% Use MATLAB's backslash operator
z_interior = A_sys \ b_sys;

% Combine with boundary conditions z_1 and z_n
z_values(1) = 0;         % Natural spline condition S''(x_1) = 0 => z_1 = 0
z_values(2:n-1) = z_interior;
z_values(n) = 0;         % Natural spline condition S''(x_n) = 0 => z_n = 0


% --- Calculate coefficients in x^3, x^2, x, 1 basis for each segment ---
% For each segment k=1, ..., num_segments, S_k(x) = A_k x^3 + B_k x^2 + C_k x + D_k
coefficients = zeros(num_segments, 4); % Columns for A_k, B_k, C_k, D_k

for k = 1:num_segments
    h_k = h_values(k);
    m_k = m_values(k); % Slope between y_k and y_{k+1}
    z_k = z_values(k);
    z_kplus1 = z_values(k+1);
    x_k = x_values(k);
    y_k = y_values(k);

    % Formulas for coefficients in (x-xk)^3, (x-xk)^2, (x-xk), 1 basis: a, b, c, d
    ak_xk = (z_kplus1 - z_k) / (6 * h_k);       % Coefficient of (x-xk)^3
    bk_xk = z_k / 2;                           % Coefficient of (x-xk)^2
    ck_xk = m_k - (h_k / 6) * (z_kplus1 + 2 * z_k); % Coefficient of (x-xk)
    dk_xk = y_k;                               % Constant term

    % Convert to x^3, x^2, x, 1 basis coefficients (A_k, B_k, C_k, D_k)
    % If S_k(x) = a(x-xk)^3 + b(x-xk)^2 + c(x-xk) + d, then Ax^3 + Bx^2 + Cx + D
    % A = a
    % B = b - 3axk
    % C = c + 3axk^2 - 2bxk
    % D = d - axk^3 + bkx^2 - cxk

    A_k = ak_xk;
    B_k = bk_xk - 3 * ak_xk * x_k;
    C_k = ck_xk + 3 * ak_xk * x_k^2 - 2 * bk_xk * x_k;
    D_k = dk_xk - ak_xk * x_k^3 + bk_xk * x_k^2 - ck_xk * x_k;


    coefficients(k, :) = [A_k, B_k, C_k, D_k]; % Store A_k, B_k, C_k, D_k

    % Construct polynomial string
    spline_strings{k} = cubic_poly_to_string(A_k, B_k, C_k, D_k); % Use helper function
end

% --- Display results in the order of the document ---

fprintf('Resultados:\n\n'); % Encabezado "Resultados:"

fprintf('Coeficientes de los trazadores:\n'); % Encabezado
% Print the coefficients matrix row by row
for row = 1:num_segments
    fprintf('%10.6f %10.6f %10.6f %10.6f\n', coefficients(row, 1), coefficients(row, 2), coefficients(row, 3), coefficients(row, 4)); % Print A_k, B_k, C_k, D_k
end

fprintf('\nTrazadores:\n'); % Encabezado
% Print each spline polynomial string
for k = 1:num_segments
    fprintf('%s\n', spline_strings{k});
end


end

% Helper function to convert polynomial coefficients [A, B, C, D] to a formatted string Ax^3 + Bx^2 + Cx + D
% A, B, C, D are coefficients of x^3, x^2, x, 1 respectively.
function poly_str = cubic_poly_to_string(A, B, C, D)
    % Build terms and concatenate with signs
    terms = {};

    % Term Ax^3
    if abs(A) >= 1e-10
        terms{end+1} = [sprintf('%.6f', A) 'x^3'];
    elseif abs(B) < 1e-10 && abs(C) < 1e-10 && abs(D) < 1e-10 % All coefficients are effectively zero
        terms{end+1} = sprintf('%.6f', 0.0);
    else % A is zero, but not all coefficients are zero, skip this term for formatting
        % Do nothing
    end

    % Term Bx^2
    if abs(B) >= 1e-10
        if B > 0 && ~isempty(terms) && ~strcmp(terms{1}, sprintf('%.6f', 0.0)) % Add '+' if B>0 and previous non-zero terms exist
             terms{end+1} = '+';
        end
        terms{end+1} = [sprintf('%.6f', B) 'x^2'];
    end

    % Term Cx
    if abs(C) >= 1e-10
        if C > 0 && ~isempty(terms) && ~strcmp(terms{1}, sprintf('%.6f', 0.0))
             terms{end+1} = '+';
        end
        terms{end+1} = [sprintf('%.6f', C) 'x'];
    end

    % Term D
    if abs(D) >= 1e-10
        if D > 0 && ~isempty(terms) && ~strcmp(terms{1}, sprintf('%.6f', 0.0))
             terms{end+1} = '+';
        end
        terms{end+1} = sprintf('%.6f', D);
    elseif isempty(terms) % This case should be caught by the all-zero check for A
         % Do nothing (should not be reached)
    end

    % Concatenate terms
    poly_str = strjoin(terms, '');

    % Handle case where the first non-zero term is positive and does not have a leading '+'.
    % strjoin handles this correctly.

end
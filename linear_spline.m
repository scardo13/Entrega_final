function [coefficients, spline_strings] = linear_spline(x_values, y_values)
% linear_spline Implements linear spline interpolation.
% Calculates and displays coefficients and polynomial strings for each linear segment.
% Matches the output format of the document for Method 12.
%
% Entradas:
%   x_values - Vector de valores de x.
%   y_values - Vector de valores de y correspondientes.
%
% Salidas:
%   coefficients   - Matrix of coefficients [a_k, b_k] for each segment.
%   spline_strings - Cell array of polynomial strings for each segment.

n = length(x_values); % Número de puntos
if length(y_values) ~= n
    error('Los vectores x y y deben tener la misma longitud.');
end

num_segments = n - 1; % Número de segmentos lineales

if num_segments <= 0
    error('Se necesitan al menos 2 puntos para crear un trazador lineal.');
end

% Inicializar arrays/celdas para almacenar coeficientes y cadenas de polinomios
coefficients = zeros(num_segments, 2); % Columna 1 para 'a_k', Columna 2 para 'b_k'
spline_strings = cell(num_segments, 1);

% Calcular coeficientes y cadenas para cada segmento
for k = 1:num_segments % k representa el segmento (desde 1 hasta num_segments)
    % Obtener puntos para el segmento actual (xk, yk) y (xk+1, yk+1)
    x_k = x_values(k);
    y_k = y_values(k);
    x_kplus1 = x_values(k+1);
    y_kplus1 = y_values(k+1);

    % Verificar si el denominador es cero (puntos x duplicados)
    if (x_kplus1 - x_k) == 0
        error('Puntos x duplicados detectados entre x(%d) y x(%d). No se puede calcular el trazador lineal.', k, k+1);
    end

    % Calcular la pendiente (a_k)
    a_k = (y_kplus1 - y_k) / (x_kplus1 - x_k);

    % Calcular el intercepto (b_k) usando el punto (xk, yk)
    b_k = y_k - a_k * x_k;

    % Almacenar coeficientes
    coefficients(k, :) = [a_k, b_k];

    % Construir la cadena del polinomio S_k(x) = a_k x + b_k
    Sk_str = '';

    % Añadir término 'ax'
    if abs(a_k) < 1e-10 % Si la pendiente es efectivamente cero
        % No añadir término 'ax', solo el término constante 'b'
    else
         Sk_str = [Sk_str sprintf('%.6f', a_k) 'x']; %#ok<AGROW>
    end

    % Añadir término 'b'
    if abs(b_k) < 1e-10 % Si el intercepto es efectivamente cero
         % Si la pendiente NO fue cero, añadir '+0.000000'
         if abs(a_k) >= 1e-10
              Sk_str = [Sk_str '+0.000000']; %#ok<AGROW>
         else % Ambos a_k y b_k son cero
              Sk_str = '0.000000'; % El polinomio es simplemente 0
         end
    else % Si el intercepto no es cero
        if b_k > 0
             Sk_str = [Sk_str '+' sprintf('%.6f', b_k)]; %#ok<AGROW>
        else % b_k < 0
             Sk_str = [Sk_str sprintf('%.6f', b_k)]; % sprintf('%.6f', b_k) incluye el signo '-', e.g., "-7.000000"
        end
    end
     % Manejar el caso donde solo el término constante existe y es positivo (asegurarse de que no empiece con + si es el único término)
     if abs(a_k) < 1e-10 && b_k > 0 && startsWith(Sk_str, '+') && length(Sk_str) > 1
          % Si solo hay un término constante positivo y empieza con +, remover el +
          Sk_str = Sk_str(2:end);
     end


    % Almacenar la cadena del polinomio
    spline_strings{k} = Sk_str;
end % Fin del bucle for k

% --- Mostrar resultados en el orden del documento ---

fprintf('Resultados:\n\n'); % Encabezado "Resultados:"

fprintf('Coeficientes de los trazadores:\n'); % Encabezado
% Imprimir la matriz de coeficientes fila por fila
for row = 1:num_segments
    fprintf('%10.6f %10.6f\n', coefficients(row, 1), coefficients(row, 2)); % Imprimir a_k y b_k
end % Fin del bucle for row (print coefficients)

fprintf('\nTrazadores:\n'); % <-- Esta es la línea 93 en mi conteo actual, no debería dar error aquí
% Imprimir cada cadena de polinomio de trazador
for k = 1:num_segments
    fprintf('%s\n', spline_strings{k});
end % Fin del bucle for k (print spline strings)

end % Fin de la función
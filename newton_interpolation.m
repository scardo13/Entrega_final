function [divided_differences_table, coefficients, poly_str] = newton_interpolation(x_values, y_values)
% newton_interpolation Calcula el polinomio interpolador usando el método de Newton (diferencias divididas).
% Muestra la tabla de diferencias divididas, coeficientes y el polinomio en el orden del documento.
%
% Entradas:
%   x_values - Vector de valores de x.
%   y_values - Vector de valores de y correspondientes.
%
% Salidas:
%   divided_differences_table - La tabla de diferencias divididas (en formato estándar).
%   coefficients              - El vector de coeficientes del polinomio de Newton.
%   poly_str                  - Una cadena de texto que representa el polinomio de Newton.

n = length(x_values);
m = length(y_values);

if n ~= m
    error('Los vectores x y y deben tener la misma longitud.');
end

% Inicializar la tabla de diferencias divididas en el formato estándar (para cálculo)
% La primera columna es la de los valores de y
ddt_standard = zeros(n, n);
ddt_standard(:, 1) = y_values(:); % Aseguramos que y_values sea columna

% Calcular las diferencias divididas (llenando la tabla estándar)
for j = 2:n % Columna (orden de la diferencia), desde la 2da columna
    for i = 1:n - j + 1 % Fila, calculamos hasta donde es posible en esa columna
         % Fórmula: f[xi...xi+k] = (f[xi+1...xi+k] - f[xi...xi+k-1]) / (xi+k - xi)
         % Indices en la matriz estándar: ddt_standard(i, j) = (ddt_standard(i+1, j-1) - ddt_standard(i, j-1)) / (x_values(i+j-1) - x_values(i))
         ddt_standard(i, j) = (ddt_standard(i + 1, j - 1) - ddt_standard(i, j - 1)) / (x_values(i + j - 1) - x_values(i));
    end
end

% --- Mostrar resultados en el orden del documento ---

fprintf('Resultados:\n\n'); % Encabezado "Resultados:"

fprintf('Tabla de diferencias divididas:\n'); % Encabezado "Tabla de diferencias divididas:"
% Imprimir la tabla en el orden del documento
% El elemento en la fila i, columna j del documento es ddt_standard(i - j + 1, j)
for i = 1:n % Fila en el documento (1 a n)
    for j = 1:n % Columna en el documento (1 a n)
        % Calcular el índice de fila correspondiente en la matriz estándar
        standard_row_index = i - j + 1;

        % Verificar si el índice de fila es válido en la matriz estándar (evitar imprimir fuera de los cálculos)
        if standard_row_index >= 1 && standard_row_index <= n
             fprintf('%10.6f ', ddt_standard(standard_row_index, j)); % Imprimir el valor correspondiente
        else
             fprintf('%10.6f ', 0.0); % Si el índice no es válido, imprimir 0 (como en la parte superior triangular del documento)
        end
    end
    fprintf('\n');
end

% Los coeficientes del polinomio de Newton son los elementos de la primera fila de la tabla ESTÁNDAR
coefficients = ddt_standard(1, :);

fprintf('\nCoeficientes del polinomio de Newton:\n'); % Encabezado "Coeficientes del polinomio de Newton:"
% Imprimir los coeficientes
for col = 1:n
    fprintf('%10.6f ', coefficients(col)); % Formato con 6 decimales
end
fprintf('\n');

% Construir la cadena del polinomio de Newton (usando los coeficientes de la tabla estándar)
poly_str = sprintf('%.6f', coefficients(1)); % El primer término es solo el primer coeficiente

for k = 2:n % Para los términos desde el segundo en adelante (k-1 es el orden del término)
    coef = coefficients(k); % El coeficiente del término de orden k-1 es ddt_standard(1, k)

     % Ignorar términos con coeficientes muy pequeños (cercanos a cero)
    if abs(coef) < 1e-10
        continue;
    end

    % Añadir el signo y el coeficiente (solo si el coeficiente no es 0)
    if coef > 0
        poly_str = [poly_str '+']; %#ok<AGROW>
    else % coef < 0
        poly_str = [poly_str '-']; %#ok<AGROW>
        coef = abs(coef); % Usamos el valor absoluto para imprimir el número después del signo
    end
     % Incluir el coeficiente numérico
     poly_str = [poly_str sprintf('%.6f', coef)]; %#ok<AGROW>


    % Añadir los términos (x - xi)
    for i = 1:k-1 % El término de orden k-1 tiene k-1 factores (x - x0)...(x - xk-2)
         % Determinar el signo y el valor de x_values(i-1) para la cadena (los puntos se usan desde x0)
         point_value = x_values(i); % i va de 1 a k-1, usamos x_values(1) a x_values(k-1)
         if point_value < 0
             % Si xi es negativo, el factor es (x - (-|xi|)) = (x + |xi|)
             poly_str = [poly_str '(x+' sprintf('%.6f', abs(point_value)) ')']; %#ok<AGROW>
         elseif point_value > 0
              % Si xi es positivo, el factor es (x - xi)
              poly_str = [poly_str '(x-' sprintf('%.6f', point_value) ')']; %#ok<AGROW> % LINEA CORREGIDA AQUÍ
         else % point_value es 0
              % Si xi es 0, el factor es (x - 0) = (x)
              poly_str = [poly_str '(x)']; %#ok<AGROW>
         end
    end
end

% Ajuste final: Si el primer término es positivo y la cadena comienza con '+', removerlo.
if startsWith(poly_str, '+')
    poly_str = poly_str(2:end);
end


fprintf('\nPolinomio de Newton:\n'); % Encabezado "Polinomio de Newton:"
fprintf('%s\n', poly_str);


end
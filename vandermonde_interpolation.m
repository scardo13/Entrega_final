function [V, a, poly_str] = vandermonde_interpolation(x_values, y_values)
% vandermonde_interpolation Calcula el polinomio interpolador usando el método de Vandermonde.
%
% Entradas:
%   x_values - Vector de valores de x.
%   y_values - Vector de valores de y correspondientes.
%
% Salidas:
%   V        - La matriz de Vandermonde construida.
%   a        - El vector de coeficientes del polinomio.
%   poly_str - Una cadena de texto que representa el polinomio interpolador.

n = length(x_values);
m = length(y_values);

if n ~= m
    error('Los vectores x y y deben tener la misma longitud.');
end

% Construir la matriz de Vandermonde
% La forma es V(i, j) = x_values(i)^(n-j) para el polinomio a_n*x^n + ... + a_1*x + a_0
% o V(i, j) = x_values(i)^(j-1) para el polinomio a_0 + a_1*x + ... + a_n*x^n
% El documento muestra la forma a3*x^3 + a2*x^2 + a1*x + a0,
% lo que implica que las columnas de la matriz son x^3, x^2, x^1, x^0.
V = zeros(n, n);
for i = 1:n
    for j = 1:n
        V(i, j) = x_values(i)^(n-j);
    end
end

% Mostrar la Matriz de Vandermonde
fprintf('Resultados:\n\n');
fprintf('Matriz de Vandermonde:\n');
for row = 1:n
    for col = 1:n
        fprintf('%10.6f ', V(row, col));
    end
    fprintf('\n');
end

% Resolver el sistema lineal Va = y
% Usamos el operador \ de MATLAB que usa métodos numéricamente estables
a = V \ y_values(:); % Aseguramos que y_values sea columna

% Mostrar los Coeficientes del polinomio
fprintf('\nCoeficientes del polinomio:\n');
for row = 1:n
    fprintf('%10.6f\n', a(row));
end

% Construir la cadena del polinomio
poly_str = '';
for j = 1:n
    coef = a(j);
    power = n - j;

    % Ignorar términos con coeficientes muy pequeños (cercanos a cero)
    if abs(coef) < 1e-10
        continue;
    end

    % Añadir el signo
    if coef > 0 && j > 1
        poly_str = [poly_str '+']; %#ok<AGROW>
    end
    if coef < 0
        poly_str = [poly_str '-']; %#ok<AGROW>
        coef = abs(coef);
    end

    % Añadir el coeficiente (si no es 1 para x^power salvo si power es 0)
    if abs(coef - 1) > 1e-10 || power == 0
         poly_str = [poly_str sprintf('%.6f', coef)]; %#ok<AGROW>
    end


    % Añadir la parte de x
    if power > 1
        poly_str = [poly_str 'x^' num2str(power)]; %#ok<AGROW>
    elseif power == 1
        % Si el coeficiente es 1 o -1, y no lo imprimimos antes, solo ponemos 'x'
        if abs(coef - 1) < 1e-10 && j > 1 && ~isempty(poly_str) && (poly_str(end) == '+' || poly_str(end) == '-')
             % Ya se agregó el signo, solo añadimos 'x' si el coeficiente es 1
             poly_str = [poly_str 'x'];%#ok<AGROW>
        elseif abs(coef - 1) > 1e-10 || j == 1 % Si el coeficiente no es 1 o es el primer término
            poly_str = [poly_str 'x'];%#ok<AGROW>
        elseif j == 1 && abs(coef - 1) < 1e-10
             poly_str = [poly_str 'x'];%#ok<AGROW>
        end

         % Si el coeficiente es -1 y no lo imprimimos antes
         if abs(coef - 1) < 1e-10 && coef > 0 && j > 1 && ~isempty(poly_str) && poly_str(end) == '-' % Si el coeficiente original era -1
             poly_str = [poly_str 'x'];%#ok<AGROW>
         end


    end
end

% Mostrar el Polinomio
fprintf('\nPolinomio:\n');
fprintf('%s\n', poly_str);


end
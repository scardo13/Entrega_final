function [x, iter, error_hist, T, C, spectral_radius] = jacobi_method(A, b, x0, tol, nmax)
% jacobi_method Implementa el método iterativo de Jacobi para resolver Ax=b.
% Calcula y muestra las matrices T, C, el radio espectral y la tabla de iteraciones.
% La matriz T calculada coincide con la definición estándar T = -inv(D)*(A-D).
%
% Entradas:
%   A     - Matriz de coeficientes.
%   b     - Vector de términos independientes.
%   x0    - Vector de aproximación inicial.
%   tol   - Tolerancia para el criterio de parada (norma 2 del error relativo o absoluto).
%   nmax  - Número máximo de iteraciones.
%
% Salidas:
%   x               - Vector solución aproximada.
%   iter            - Número de iteraciones realizadas.
%   error_hist      - Historial de errores en cada iteración.
%   T               - Matriz de iteración T del método Jacobi (-inv(D)*(A-D)).
%   C               - Vector C del método Jacobi (inv(D)*b).
%   spectral_radius - Radio espectral de la matriz T.

[n, m] = size(A);
if n ~= m
    error('La matriz A debe ser cuadrada.');
end
if size(b, 1) ~= n || size(b, 2) ~= 1
    error('El vector b debe ser un vector columna del mismo tamaño que A.');
end
if size(x0, 1) ~= n || size(x0, 2) ~= 1
    error('El vector inicial x0 debe ser un vector columna del mismo tamaño que A.');
end

% --- Cálculo de T, C y Radio Espectral ---
% A = D + L_strict + U_strict
D = diag(diag(A));        % Matriz diagonal de A

% Verificar si algún elemento diagonal es cero (necesario para D^-1)
if any(diag(D) == 0)
    error('La matriz diagonal D tiene elementos cero. El método Jacobi no es aplicable.');
end

% Calcular inv(D)
invD = diag(1./diag(D));

% Matriz de iteración T = -inv(D) * (A - D) -- Fórmula estándar
T = -invD * (A - D);

% Vector C = inv(D) * b -- Fórmula estándar
C = invD * b;

% Calcular el radio espectral de T
eigenvalues = eig(T);
spectral_radius = max(abs(eigenvalues));

% --- Mostrar T, C y Radio Espectral ---
fprintf('Resultados:\n\n'); % Encabezado "Resultados:"

fprintf('T:\n'); % Encabezado "T:"
for row = 1:n
    for col = 1:n
        fprintf('%10.6f ', T(row, col)); % Formato con 6 decimales
    end
    fprintf('\n');
end

fprintf('\nC:\n'); % Encabezado "C:"
for row = 1:n
    fprintf('%10.6f\n', C(row)); % Formato con 6 decimales
end

fprintf('\nradio espectral:\n'); % Encabezado "radio espectral:"
fprintf('%10.6f\n', spectral_radius); % Formato con 6 decimales


% --- Proceso Iterativo de Jacobi ---
x = x0;
iter = 0;
error_hist = [];

% Preparar la tabla de resultados
fprintf('\n| iter |         E         | ');
for i = 1:n
    fprintf('      x%d       | ', i);
end
fprintf('\n');
fprintf('------|-------------------|');
for i = 1:n
    fprintf('-----------------|');
end
fprintf('\n');


% Mostrar la aproximación inicial (iteración 0)
fprintf('%5d | %17.6e | ', iter, 0.0); % Error inicial se considera 0
for i = 1:n
    fprintf('%16.6f | ', x(i));
end
fprintf('\n');

% Bucle iterativo
for k = 1:nmax
    iter = k; % Usamos k como el número de iteración actual
    x_old = x; % Guardar el valor anterior para calcular el error
    x_new = zeros(n, 1); % Vector para la nueva aproximación

    % Iteración elemento a elemento
    for i = 1:n
        sum_off_diagonal = 0;
        for j = 1:n
            if i ~= j
                sum_off_diagonal = sum_off_diagonal + A(i, j) * x_old(j); % Usa valores de la iteración anterior (k)
            end
        end

        % Asegurarse de que A(i,i) no sea cero
         if A(i,i) == 0
             error('Elemento diagonal A(%d,%d) es cero. El método Jacobi no es aplicable directamente.', i, i);
         end

        x_new(i) = (b(i) - sum_off_diagonal) / A(i, i);
    end

    x = x_new; % Actualizar la aproximación

    % Calcular el error (norma 2 de la diferencia entre iteraciones)
    error = norm(x - x_old, 2);
    error_hist = [error_hist, error]; %#ok<AGROW> % Almacenar historial de errores

    % Mostrar los resultados de la iteración actual
    fprintf('%5d | %17.6e | ', iter, error);
    for i = 1:n
        fprintf('%16.6f | ', x(i));
    end
    fprintf('\n');

    % Criterio de parada
    if error < tol
        fprintf('\nConvergencia alcanzada en %d iteraciones.\n', iter);
        break;
    end

    % Si se alcanza el número máximo de iteraciones sin converger
    if k == nmax
        warning('El método no convergió después de %d iteraciones.', nmax);
    end
end

% Mostrar el resultado final del vector x
fprintf('\nDespués de aplicar sustitución progresiva y regresiva (Resultado Final Iterativo):\n'); % Ajuste en el texto para reflejar que es resultado iterativo
fprintf('X:\n');
for row = 1:n
    fprintf('%10.6f\n', x(row));
end


end
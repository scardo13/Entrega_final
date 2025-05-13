function [x, iter, error_hist, T, C, spectral_radius] = sor_method(A, b, x0, tol, omega, nmax)
% sor_method Implementa el método de sobre-relajación sucesiva (SOR) para resolver Ax=b.
% Calcula y muestra las matrices T, C y el radio espectral, además de la tabla de iteraciones.
%
% Entradas:
%   A     - Matriz de coeficientes.
%   b     - Vector de términos independientes.
%   x0    - Vector de aproximación inicial.
%   tol   - Tolerancia para el criterio de parada (norma 2 del error relativo o absoluto).
%   omega - Factor de relajación (0 < omega < 2).
%   nmax  - Número máximo de iteraciones.
%
% Salidas:
%   x               - Vector solución aproximada.
%   iter            - Número de iteraciones realizadas.
%   error_hist      - Historial de errores en cada iteración.
%   T               - Matriz de iteración T del método SOR.
%   C               - Vector C del método SOR.
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
% No detenemos la ejecución por el valor de omega, solo advertimos si está fuera del rango común.
if omega <= 0 || omega >= 2
    warning('El factor de relajación (omega) %.4f está fuera del rango (0, 2). La convergencia no está garantizada.', omega);
end

% --- Cálculo de T, C y Radio Espectral ---
% A = D - L - U
D = diag(diag(A));        % Matriz diagonal de A
L_matrix = -tril(A, -1);  % Parte estrictamente triangular inferior de A (con signo negativo)
U_matrix = -triu(A, 1);   % Parte estrictamente triangular superior de A (con signo negativo)

% Matriz M = D - omega*L
M = D - omega * L_matrix;

% Matriz N = (1-omega)*D + omega*U
N = (1 - omega) * D + omega * U_matrix;

% Verificar si M es singular antes de invertir
if det(M) == 0
    error('La matriz (D - omega*L) es singular. No se pueden calcular T y C.');
end

% Matriz de iteración T = inv(M) * N
T = inv(M) * N;

% Vector C = inv(M) * (omega*b)
C = inv(M) * (omega * b);

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


% --- Proceso Iterativo de SOR ---
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

    % Iteración elemento a elemento de SOR
    for i = 1:n
        sum_before = 0;
        for j = 1:i-1
            sum_before = sum_before + A(i, j) * x(j); % Usa valores actualizados (iteracion actual)
        end

        sum_after = 0;
        for j = i+1:n
            sum_after = sum_after + A(i, j) * x_old(j); % Usa valores de la iteración anterior
        end

        % Asegurarse de que A(i,i) no sea cero
        if A(i,i) == 0
             error('Elemento diagonal A(%d,%d) es cero. El método SOR no es aplicable directamente.', i, i);
        end

        x(i) = (1 - omega) * x_old(i) + (omega / A(i, i)) * (b(i) - sum_before - sum_after);
    end

    % Calcular el error (norma 2 de la diferencia entre iteraciones)
    error = norm(x - x_old, 2);
    error_hist = [error_hist, error]; %#ok<AGROW> % Almacenar historial de errores

    % Mostrar los resultados de la iteración actual
    fprintf('%5d | %17.6e | ', iter, error); % Usamos iter (que es k)
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
fprintf('\nDespués de aplicar sustitución progresiva y regresiva (Resultado Final Iterativo):\n'); % Ajuste en el texto para refle
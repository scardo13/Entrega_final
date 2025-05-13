function [x, iter, error_hist, T, C, spectral_radius] = gauss_seidel_method(A, b, x0, tol, nmax)
% gauss_seidel_method Implementa el método iterativo de Gauss-Seidel para resolver Ax=b.
% Calcula y muestra las matrices T, C, el radio espectral y la tabla de iteraciones.
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
%   T               - Matriz de iteración T del método Gauss-Seidel ((D-L)^-1 * U).
%   C               - Vector C del método Gauss-Seidel ((D-L)^-1 * b).
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

% Inicializar variables de salida para asegurar que estén asignadas en todas las rutas (silencia posibles advertencias/errores)
x = [];              % Inicializa x
iter = 0;            % Inicializa iter
error_hist = [];     % Inicializa error_hist
T = [];              % Inicializa T
C = [];              % Inicializa C
spectral_radius = NaN; % Inicializa spectral_radius

% --- Cálculo de T, C y Radio Espectral ---
% A = D - L - U (donde L es estrictamente triangular inferior, U es estrictamente triangular superior)
D = diag(diag(A));        % Matriz diagonal de A
L_matrix = -tril(A, -1);  % -L (parte estrictamente triangular inferior de A) -> L = -(-L)
U_matrix = -triu(A, 1);   % -U (parte estrictamente triangular superior de A) -> U = -(-U)


% M = D - L = D + tril(A, -1)
M = D + tril(A, -1); % O equivalentemente: M = tril(A);

% Verificar si M es singular antes de calcular la inversa o usar \
if det(M) == 0
    error('La matriz (D - L) es singular. El método Gauss-Seidel no es aplicable directamente.');
end

% Matriz de iteración T = (D - L)^-1 * U
T = M \ U_matrix; % Usando el operador \ para resolver M*T = U

% Vector C = (D - L)^-1 * b
C = M \ b; % Usando el operador \ para resolver M*C = b

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


% --- Proceso Iterativo de Gauss-Seidel ---
x = x0; % Asigna el valor inicial a x para el proceso iterativo
iter = 0; % Reinicia iter a 0 para el inicio de la tabla de iteraciones
error_hist = []; % Limpia error_hist para el inicio de la tabla de iteraciones


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

    % Iteración elemento a elemento de Gauss-Seidel
    for i = 1:n
        sum_before = 0;
        for j = 1:i-1
            sum_before = sum_before + A(i, j) * x(j); % Usa valores actualizados (iteracion actual)
        end

        sum_after = 0;
        for j = i+1:n
            sum_after = sum_after + A(i, j) * x_old(j); % Usa valores de la iteración anterior
        end

        % Asegurarse de que A(i,i) no sea cero (necesario para la división)
         if A(i,i) == 0
             error('Elemento diagonal A(%d,%d) es cero. El método Gauss-Seidel no es aplicable directamente.', i, i);
         end

        x(i) = (b(i) - sum_before - sum_after) / A(i, i);
    end

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
fprintf('\nDespués de aplicar el método iterativo de Gauss-Seidel:\n'); % Texto ajustado
fprintf('X:\n');
for row = 1:n
    fprintf('%10.6f\n', x(row));
end


end
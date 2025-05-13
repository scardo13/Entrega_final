function [L, U, P, x] = lu_partial_pivot(A, b)
% lu_partial_pivot Implementa la descomposición LU con pivoteo parcial para resolver Ax=b.
% Muestra el estado de la matriz y P, L, U en cada etapa.
%
% Entradas:
%   A - Matriz cuadrada de coeficientes.
%   b - Vector columna de términos independientes.
%
% Salidas:
%   L - Matriz triangular inferior.
%   U - Matriz triangular superior.
%   P - Matriz de permutación tal que PA = LU.
%   x - Vector solución del sistema Ax=b.

[n, m] = size(A);
if n ~= m
    error('La matriz A debe ser cuadrada.');
end
if size(b, 1) ~= n || size(b, 2) ~= 1
    error('El vector b debe ser un vector columna del mismo tamaño que A.');
end

% Inicializar L como matriz identidad, U como copia de A y P como matriz identidad
L = eye(n);
U = A;
P = eye(n);

% Mostrar Matriz inicial (Etapa 0) para validación
fprintf('Etapa 0:\n');
% Para mostrar la matriz con formato fijo, iteramos sobre sus elementos
for row = 1:n
    for col = 1:n
        fprintf('%10.6f ', U(row, col)); % Mostrar el estado inicial (matriz A)
    end
    fprintf('\n');
end


% Proceso de eliminación Gaussiana con pivoteo parcial
for k = 1:n-1
    % Encontrar el pivote más grande en la columna k, desde la fila k hacia abajo
    [~, pivot_row] = max(abs(U(k:n, k)));
    pivot_row = pivot_row + k - 1; % Ajustar el índice a la matriz completa

    % Si el pivote no está en la fila actual, intercambiar filas
    if pivot_row ~= k
        % Intercambiar filas en U
        U([k, pivot_row], :) = U([pivot_row, k], :);
        % Intercambiar filas en L (solo las columnas ya procesadas)
        L([k, pivot_row], 1:k-1) = L([pivot_row, k], 1:k-1);
        % Intercambiar filas en P
        P([k, pivot_row], :) = P([pivot_row, k], :);
    end

    % Verificar si el pivote es cero después del pivoteo
    if U(k, k) == 0
        error('El elemento pivote es cero incluso con pivoteo parcial. La matriz es singular.');
    end

    % Realizar la eliminación
    for i = k+1:n
        % Calcular el multiplicador
        factor = U(i, k) / U(k, k);
        L(i, k) = factor; % Almacenar el multiplicador en L

        % Eliminar el elemento debajo del pivote en U
        U(i, k:n) = U(i, k:n) - factor * U(k, k:n);
    end

    % Mostrar el estado de la matriz (U actual), L, U y P en esta etapa
    fprintf('Etapa %d:\n', k);
    % Mostrar el estado actual de U después de la eliminación de la columna k
    for row = 1:n
        for col = 1:n
            fprintf('%10.6f ', U(row, col));
        end
        fprintf('\n');
    end

    fprintf('L:\n');
     for row = 1:n
        for col = 1:n
            fprintf('%10.6f ', L(row, col));
        end
        fprintf('\n');
    end
    fprintf('U:\n');
     for row = 1:n
        for col = 1:n
            fprintf('%10.6f ', U(row, col));
        end
        fprintf('\n');
    end
     fprintf('P:\n');
     for row = 1:n
        for col = 1:n
            fprintf('%10.6f ', P(row, col));
        end
        fprintf('\n');
    end
end

% Resolver Ly = Pb usando sustitución progresiva
Pb = P * b; % Aplicar la misma permutación a b
y = zeros(n, 1);
for i = 1:n
    sum_ly = 0;
    for j = 1:i-1
        sum_ly = sum_ly + L(i, j) * y(j);
    end
    y(i) = (Pb(i) - sum_ly) / L(i, i);
end

% Resolver Ux = y usando sustitución regresiva
x = zeros(n, 1);
for i = n:-1:1
    sum_ux = 0;
    for j = i+1:n
        sum_ux = sum_ux + U(i, j) * x(j);
    end
    x(i) = (y(i) - sum_ux) / U(i, i);
end

% Mostrar el resultado final del vector x
fprintf('\nDespués de aplicar sustitución progresiva y regresiva:\n');
fprintf('X:\n');
for row = 1:n
    fprintf('%10.6f\n', x(row));
end

end
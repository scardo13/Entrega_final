function [L, U, x] = lu_simple_gauss(A, b)
% lu_simple_gauss Implementa la descomposición LU sin pivoteo para resolver Ax=b.
% Muestra el estado de la matriz en cada etapa de la eliminación.
%
% Entradas:
%   A - Matriz cuadrada de coeficientes.
%   b - Vector columna de términos independientes.
%
% Salidas:
%   L - Matriz triangular inferior.
%   U - Matriz triangular superior.
%   x - Vector solución del sistema Ax=b.

[n, m] = size(A);
if n ~= m
    error('La matriz A debe ser cuadrada.');
end
if size(b, 1) ~= n || size(b, 2) ~= 1
    error('El vector b debe ser un vector columna del mismo tamaño que A.');
end

% Inicializar L como la matriz identidad y U como una copia de A
L = eye(n);
U = A;

% Mostrar Matriz inicial (Etapa 0) para validación
fprintf('Etapa 0:\n');
% Para mostrar la matriz con formato fijo, iteramos sobre sus elementos
for row = 1:n
    for col = 1:n
        fprintf('%10.6f ', U(row, col)); % Mostrar el estado inicial (matriz A)
    end
    fprintf('\n');
end


% Proceso de eliminación Gaussiana
for k = 1:n-1
    % Verificar si el pivote es cero
    if U(k, k) == 0
        error('El elemento pivote es cero. No se puede proceder sin pivoteo.');
    end

    for i = k+1:n
        % Calcular el multiplicador
        factor = U(i, k) / U(k, k);
        L(i, k) = factor; % Almacenar el multiplicador en L

        % Eliminar el elemento debajo del pivote en U
        U(i, k:n) = U(i, k:n) - factor * U(k, k:n);
    end

    % Mostrar el estado de la matriz (U actual), L y U en esta etapa
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
end

% Resolver Ly = b usando sustitución progresiva
y = zeros(n, 1);
for i = 1:n
    sum_ly = 0;
    for j = 1:i-1
        sum_ly = sum_ly + L(i, j) * y(j);
    end
    y(i) = (b(i) - sum_ly) / L(i, i);
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
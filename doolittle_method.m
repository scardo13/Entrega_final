function [L, U, x] = doolittle_method(A, b)
% doolittle_method Implements Doolittle decomposition (A=LU, L unit lower triangular) to solve Ax=b.
% L has ones on the diagonal.
% Displays L and U at each stage and the final solution.
% Matches the output format of the document for Method 4.
%
% Entradas:
%   A - Matriz cuadrada de coeficientes.
%   b - Vector columna de términos independientes.
%
% Salidas:
%   L - Matriz triangular inferior con unos en la diagonal.
%   U - Matriz triangular superior.
%   x - Vector solución del sistema Ax=b.

[n, m] = size(A);
if n ~= m
    error('La matriz A debe ser cuadrada.');
end
if size(b, 1) ~= n || size(b, 2) ~= 1
    error('El vector b debe ser un vector columna del mismo tamaño que A.');
end

% Inicializar L como matriz identidad n x n y U como matriz n x n con ceros
% L tendrá unos en la diagonal y elementos por debajo. U tendrá elementos en
% la diagonal y por encima.
L = eye(n);
U = zeros(n);

% Mostrar Matriz inicial (Etapa 0) para validación
fprintf('Etapa 0:\n');
% Para mostrar la matriz con formato fijo, iteramos sobre sus elementos
for row = 1:n
    for col = 1:n
        fprintf('%10.6f ', A(row, col)); % Mostrar la matriz A inicial
    end
    fprintf('\n');
end

% Proceso de Doolittle decomposition (A = LU) - calculando fila i de U y columna i de L
for i = 1:n % i representa la fila de U y columna de L que se está calculando (desde 1 hasta n)
    % Compute elements of U in row i (u_ij for j >= i)
    for j = i:n % j representa la columna (desde i hasta n)
        sum_u = 0;
        for p = 1:i-1 % p representa el índice de la suma (desde 1 hasta i-1)
            sum_u = sum_u + L(i, p) * U(p, j);
        end
        U(i, j) = A(i, j) - sum_u;
    end

    % Check for zero on the diagonal of U (pivot element u_ii)
    if U(i, i) == 0
        error('Elemento diagonal U(%d,%d) es cero. La descomposición de Doolittle no es posible sin pivoteo.', i, i);
    end

    % Compute elements of L in column i (l_ik for k > i)
    % L(i, i) is already 1 due to initialization
    for k = i+1:n % k representa la fila (desde i+1 hasta n)
        sum_l = 0;
        for p = 1:i-1 % p representa el índice de la suma (desde 1 hasta i-1)
            sum_l = sum_l + L(k, p) * U(p, i);
        end
        L(k, i) = (A(k, i) - sum_l) / U(i, i);
    end

    % Mostrar matrices L y U en esta etapa (después de calcular fila i de U y columna i de L)
    fprintf('Etapa %d:\n', i);
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
    % L tiene unos en la diagonal, por lo tanto L(i, i) es 1, no se necesita división por L(i,i)
    y(i) = (b(i) - sum_ly); % Equivalente a (b(i) - sum_ly) / L(i,i)
end

% Resolver Ux = y usando sustitución regresiva
x = zeros(n, 1);
for i = n:-1:1
    sum_ux = 0;
    for j = i+1:n
        sum_ux = sum_ux + U(i, j) * x(j);
    end
    % U(i, i) es un elemento diagonal de U, el cual no debe ser cero (verificado durante la descomposición)
     if U(i, i) == 0
         error('Error interno: Elemento diagonal U(%d,%d) es cero durante la sustitución regresiva.', i, i);
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
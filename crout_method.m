function [L, U, x] = crout_method(A, b)
% crout_method Implements Crout decomposition to solve Ax=b.
% U has ones on the diagonal.
% Displays L and U at each stage and the final solution.
% Matches the output format of the document for Method 3.
%
% Entradas:
%   A - Matriz cuadrada de coeficientes.
%   b - Vector columna de términos independientes.
%
% Salidas:
%   L - Matriz triangular inferior.
%   U - Matriz triangular superior con unos en la diagonal.
%   x - Vector solución del sistema Ax=b.

[n, m] = size(A);
if n ~= m
    error('La matriz A debe ser cuadrada.');
end
if size(b, 1) ~= n || size(b, 2) ~= 1
    error('El vector b debe ser un vector columna del mismo tamaño que A.');
end

% Inicializar L como matriz n x n con ceros y U como matriz identidad n x n
% L tendrá elementos en la diagonal y por debajo. U tendrá unos en la diagonal
% y elementos por encima.
L = zeros(n);
U = eye(n);

% Mostrar Matriz inicial (Etapa 0) para validación
fprintf('Etapa 0:\n');
% Para mostrar la matriz con formato fijo, iteramos sobre sus elementos
for row = 1:n
    for col = 1:n
        fprintf('%10.6f ', A(row, col)); % Mostrar la matriz A inicial
    end
    fprintf('\n');
end

% Proceso de Crout decomposition
for k = 1:n % k representa la columna de L y fila de U que se está calculando
    % Compute elements of L in column k (l_ik for i >= k)
    for i = k:n
        sum_l = 0;
        for p = 1:k-1
            sum_l = sum_l + L(i, p) * U(p, k);
        end
        L(i, k) = A(i, k) - sum_l;
    end

    % Check for zero on the diagonal of L (pivot element)
    if L(k, k) == 0
        error('Elemento diagonal L(%d,%d) es cero. La descomposición de Crout no es posible sin pivoteo.', k, k);
    end

    % Compute elements of U in row k (u_kj for j > k)
    % U(k, k) is already 1 due to initialization
    for j = k+1:n
        sum_u = 0;
        for p = 1:k-1
            sum_u = sum_u + L(k, p) * U(p, j);
        end
        U(k, j) = (A(k, j) - sum_u) / L(k, k);
    end

    % Mostrar matrices L y U en esta etapa
    fprintf('Etapa %d:\n', k);
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
    % L(i, i) is a diagonal element of L, which should not be zero
    if L(i, i) == 0
         error('Error interno: Elemento diagonal L(%d,%d) es cero durante la sustitución progresiva.', i, i);
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
    % U has ones on the diagonal, so U(i, i) is 1, no division by U(i,i) is needed
    x(i) = (y(i) - sum_ux); % Divided by U(i,i) which is 1
end

% Mostrar el resultado final del vector x
fprintf('\nDespués de aplicar sustitución progresiva y regresiva:\n');
fprintf('X:\n');
for row = 1:n
    fprintf('%10.6f\n', x(row));
end


end
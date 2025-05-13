function [L, U, x] = cholesky_method(A, b)
% cholesky_method Implements Cholesky decomposition (A=LL^T) to solve Ax=b.
% Displays L and U=L^T at each stage and the final solution.
% Matches the output format of the document for Method 5, including complex numbers.
%
% Entradas:
%   A - Matriz cuadrada de coeficientes.
%   b - Vector columna de términos independientes.
%
% Salidas:
%   L - Matriz triangular inferior.
%   U - Matriz triangular superior (L^T).
%   x - Vector solución del sistema Ax=b.

[n, m] = size(A);
if n ~= m
    error('La matriz A debe ser cuadrada.');
end
if size(b, 1) ~= n || size(b, 2) ~= 1
    error('El vector b debe ser un vector columna del mismo tamaño que A.');
end

% Inicializar L como una matriz n x n de ceros
L = zeros(n);

% Mostrar Matriz inicial (Etapa 0) para validación
fprintf('Etapa 0:\n');
% Usar impresión personalizada para la matriz inicial para manejar números complejos y formato
for row = 1:n
    for col = 1:n
        val = A(row, col);
        if isreal(val) % Si es real, imprimir con .6f
            fprintf('%10.6f ', val);
        else % Si es complejo
            real_part = real(val);
            imag_part = imag(val);
            if real_part == 0 % Si la parte real es cero, imprimir solo la parte imaginaria con 'i'
                fprintf('%10.6fi ', imag_part);
            elseif imag_part == 0 % Si la parte imaginaria es cero (aunque isreal ya lo verifica)
                 fprintf('%10.6f ', real_part);
            else % Si ambas partes son distintas de cero
                fprintf('%10.6f%+.6fi ', real_part, imag_part);
            end
        end
    end
    fprintf('\n');
end


% Proceso de Cholesky decomposition (A = LL^T) - calculando L columna por columna
for i = 1:n % i representa la columna de L que se está calculando (desde 1 hasta n)
    % Calcular el elemento diagonal l_ii
    sum_diag = 0;
    for p = 1:i-1
        % En la fórmula A=LL^T, el término es la suma de los cuadrados de los elementos l_ip
        % de la fila i en las columnas anteriores a i.
        sum_diag = sum_diag + L(i, p) * L(i, p);
    end
    % El término dentro de la raíz cuadrada es a_ii - sum_diag
    val_inside_sqrt = A(i, i) - sum_diag;

    % Calcular l_ii como la raíz cuadrada del término. MATLAB's sqrt() maneja números complejos.
    L(i, i) = sqrt(val_inside_sqrt);

    % Verificar si el elemento diagonal L(i,i) es cero (con una pequeña tolerancia numérica)
    if abs(L(i, i)) < 1e-12
        error('Elemento diagonal L(%d,%d) es cero. La descomposición de Cholesky no es posible (matriz no definida positiva o singular).', i, i);
    end

    % Calcular los elementos fuera de la diagonal l_ki para k > i (en la columna i)
    for k = i+1:n % k representa la fila (desde i+1 hasta n)
        sum_off_diag = 0;
        for p = 1:i-1
             % En la fórmula A=LL^T, a_ki = sum_{p=1}^{i-1} l_kp * l_ip + l_ki * l_ii
             % Reorganizando para l_ki: l_ki = (a_ki - sum_{p=1}^{i-1} l_kp * l_ip) / l_ii
             sum_off_diag = sum_off_diag + L(k, p) * L(i, p);
        end
         L(k, i) = (A(k, i) - sum_off_diag) / L(i, i);
    end

    % Establecer U = L.' (L transpuesta no conjugada) después de calcular la columna i de L
    % Esto coincide con la forma en que U se presenta en el documento como la transpuesta de L.
    U = L.';


    % Mostrar matrices L y U=L^T en esta etapa (después de calcular la columna i de L)
    fprintf('Etapa %d:\n', i);
    fprintf('L:\n');
    % Impresión personalizada para matrices con números complejos
    for row = 1:n
        for col = 1:n
            val = L(row, col);
            if isreal(val)
                fprintf('%10.6f ', val);
            else
                real_part = real(val);
                imag_part = imag(val);
                if real_part == 0 && imag_part ~= 0
                     fprintf('%10.6fi ', imag_part); % Formato como "imag i"
                elseif imag_part == 0 && real_part ~= 0
                     fprintf('%10.6f ', real_part); % Formato como "real" (aunque isreal ya lo verifica)
                else
                    fprintf('%10.6f%+.6fi ', real_part, imag_part); % Formato como "real + imag i"
                end
            end
        end
        fprintf('\n');
    end

    fprintf('U:\n');
     for row = 1:n
        for col = 1:n
            val = U(row, col);
             if isreal(val)
                fprintf('%10.6f ', val);
            else
                real_part = real(val);
                imag_part = imag(val);
                if real_part == 0 && imag_part ~= 0
                     fprintf('%10.6fi ', imag_part);
                elseif imag_part == 0 && real_part ~= 0
                     fprintf('%10.6f ', real_part);
                else
                    fprintf('%10.6f%+.6fi ', real_part, imag_part);
                end
            end
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
    % L(i, i) es un elemento diagonal de L, no debe ser cero si la matriz no es singular
     if abs(L(i, i)) < 1e-12
         error('Error interno: Elemento diagonal L(%d,%d) es cero durante la sustitución progresiva.', i, i);
     end
    y(i) = (b(i) - sum_ly) / L(i, i);
end

% Resolver U x = y usando sustitución regresiva (U = L^T)
x = zeros(n, 1);
for i = n:-1:1
    sum_ux = 0;
    for j = i+1:n
        sum_ux = sum_ux + U(i, j) * x(j);
    end
    % U(i, i) es el elemento diagonal de U (que es L^T), por lo tanto U(i,i) = L(i,i)
     if abs(U(i, i)) < 1e-12
         error('Error interno: Elemento diagonal U(%d,%d) es cero durante la sustitución regresiva.', i, i);
     end
    x(i) = (y(i) - sum_ux) / U(i, i); % Dividido por U(i,i) que es L(i,i)
end

% Mostrar el resultado final del vector x
fprintf('\nDespués de aplicar sustitución progresiva y regresiva:\n');
fprintf('X:\n');
for row = 1:n
    fprintf('%10.6f\n', x(row));
end


end
function [Q,R] = gram_schmidt(A)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);
    for j = 1:n
        v = A(:, j);
        for i = 1:(j - 1)
            R(i, j) = A(:, j)' * Q(:, i);
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end
end

function [Q,R] = gram_schmidt_modificado(A)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);
    for j = 1:n
        v = A(:, j);
        for i = 1:(j - 1)
            R(i, j) = v' * Q(:, i);
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end
endfunction

function [Q,R] = gram_schmidt_modificado(A)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n, n);
    for j = 1:n
        v = A(:, j);
        for i = 1:(j - 1)
            R(i, j) = v' * Q(:, i);
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end
endfunction

function [U, R] = householder(A)
    [m, n] = size(A)
    U = zeros(m, n)
    R = A
    for i = 1:(n - 1)
        v = R(i:$, i)
        v(1) = v(1) - norm(v) * sign(v(1))
        u = v / norm(v)
        U(i:$, i) = u
        H = eye(size(u)(1),size(u)(1)) - 2*u*u'
        Q = eye(m, m)
        Q(i:$, i:$) = H
        R = Q * R
    end
endfunction

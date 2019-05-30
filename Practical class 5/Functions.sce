function [Q,R] = gram_schmidt(A)
    [m, n] = size(A);
    //Initialization of Q and R matrices
    Q = zeros(m, n);
    R = zeros(n, n);
    //For loop to find Q and R matrices
    for j = 1:n
        v = A(:, j);
        for i = 1:(j - 1)
            R(i, j) = A(:, j)' * Q(:, i);
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        //Q columns
        Q(:, j) = v / R(j, j);
    end
end

function [Q,R] = gram_schmidt_modificado(A)
    [m, n] = size(A);
    //Initialization of Q and R matrices
    Q = zeros(m, n);
    R = zeros(n, n);
    //For loop to find Q and R matrices
    for j = 1:n
        v = A(:, j);
        for i = 1:(j - 1)
            //The modification:
            R(i, j) = v' * Q(:, i);
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v);
        //Q columns
        Q(:, j) = v / R(j, j);
    end
endfunction

function [U, R] = householder(A)
    [m, n] = size(A)
    //Initialization of U and R matrices
    U = zeros(m, n)
    R = A
    for i = 1:n
        v = R(i:$, i)
        //Here we choose the fasthest Hx from a_i depending on its signal
        v(1) = v(1) - norm(v) * sign(v(1))
        //If v ~= 0
        if norm(v) ~= 0
        //Unit vector
            u = v / norm(v)
        else
            u = v
        end
        //Save u in U
        U(i:$, i) = u
        //Find H  
        H = eye(size(u)(1),size(u)(1)) - 2*u*u'
        //Start finding Q
        Q = eye(m, m)
        //Find Q
        Q(i:$, i:$) = H
        //Update R matrix
        R = Q * R
    end
    //What should be zero will be zero
    for i = 1:n
        R((i + 1):$, i) = zeros(m - i, 1)
    end
endfunction

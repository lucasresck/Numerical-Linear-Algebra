function [xk, normdif, k, residue] = jacobi(A, b, x0, E, M, norma)
    dim = size(A)(1);
    //Obtém L + U
    LU = zeros(dim, dim);
    for i = 1:dim
        for j = 1:(i - 1)
            LU(i, j) = A(i, j);
        end
        for j = (i + 1):dim
            LU(i, j) = A(i, j);
        end
    end
    //Obtém inv(D)
    invD = zeros(dim, dim);
    for i = 1:dim
        invD(i, i) = 1/A(i, i);
    end
    //Calcula as matrizes do método
    Mj = - invD * LU;
    cj = invD * b;
    k = 1;
    xant = x0;
    normdif = E;
    //Itera
    while (k <= M && normdif >= E)
        x = Mj * xant + cj;
        //Calcula a norma da diferença
        normdif = norm(x - xant, norma);
        xant = x;
        k = k + 1;
    end
    if k > M then
        disp("k > M");
    else
        disp("||xk - xk-1|| < E");
    end
    k = k - 1;
    xk = x;
    //Calcula a norma do resíduo
    residue = norm(b - A * xk, norma)
endfunction

function [xk, normdif, k, residue] = gauss_seidel(A, b, x0, E, M, norma)
    dim = size(A)(1);
    //Obtém L + D
    LD = zeros(dim, dim);
    for i = 1:dim
        for j = 1:i
            LD(i, j) = A(i, j);
        end
    end
    //Obtém U
    U = zeros(dim, dim);
    for i = 1:dim
        for j = (i + 1):dim
            U(i, j) = A(i, j);
        end
    end
    //Calcula as matrizes do método
    Mg = - inv(LD) * U;
    cg = inv(LD) * b;
    k = 1;
    xant = x0;
    normdif = E;
    //Itera
    while (k <= M && normdif >= E)
        x = Mg * xant + cg;
        //Calcula a norma da diferença
        normdif = norm(x - xant, norma);
        xant = x;
        k = k + 1;
    end
    if k > M then
        disp("k > M");
    else
        disp("||xk - xk-1|| < E");
    end
    k = k - 1;
    xk = x;
    //Calcula a norma do resíduo
    residue = norm(b - A * xk, norma)
endfunction

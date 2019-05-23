function [lambda, x, N] = Metodo_potencia(A, x0, epsilon, M)
    //Método da potência para cálculo de aproximações do módulo do autovalor de maior módulo e do autovetor associado a ele.
    p = find(abs(x0) == norm(x0, %inf))(1,1);   //Encontra-se o menor índice tal que sua coordenada tenha módulo igual à norma infinito
    x = x0 / x0(p, 1); //Divide-se pela coordenada de módulo igual à norma infinito
    N = 1;
    while N <= M    //Até exceder o número máximo de iterações
        y = A * x;  //Iteração
        lambda = y(p, 1);
        if norm(y, %inf) == 0 then  //Se a norma infinito de y for 0, então A tem autovalor 0
            disp("A tem autovalor 0.");
            break;
        end
        y = y / y(p, 1);
        err = norm(x - y, %inf);    //Erro
        x = y;
        if err < epsilon    //Para a iteração se o erro for menor do que a tolerância
            disp("Procedimento bem-sucedido.");
            break;
        end
        N = N + 1;
    end
    if N > M then
        N = N - 1;
        disp("Número máximo de iterações excedido.");
    end
endfunction

function [lambda, x, N] = Metodo_potencia_simetrico(A, x0, epsilon, M)
    //Método da potência simétrico para cálculo de aproximações do módulo do autovalor de maior módulo e do autovetor associado a ele para matrizes simétricas.
    x = x0 / norm(x0, 2);
    N = 1;
    while N <= M
        y = A * x;
        lambda = x'*y;
        if norm(y, 2) == 0 then
            disp("A tem autovalor 0.");
            break;
        end
        y = y / norm(y, 2);
        err1 = norm(x - y, 2);
        err2 = norm(x + y, 2);
        x = y;
        if err1 < epsilon || err2 < epsilon then
            disp("Procedimento bem-sucedido.");
            break;
        end
        N = N + 1;
    end
    if N > M then
        N = N - 1;
        disp("Número máximo de iterações excedido.");
    end
endfunction

function [lambda, x, N] = Metodo_potencia_inversa(A, x0, epsilon, alfa, M)
    //Método da potência inversa para cálculo de aproximações do módulo do autovalor mais próximo de alfa e o autovetor relativo a ele.
    if alfa == %inf then    //Se o alfa for infinito, então calculamos uma aproximação
        alfa = x0' * A * x0 / (x' * x);
    end
    p = find(abs(x0) == norm(x0, %inf))(1,1);
    x = x0 / x0(p, 1);
    N = 1;
    while N <= M
        I = eye(size(A)(1, 1), size(A)(1, 1)); //Matriz identidade de orden igual à de A
        y = linsolve(A - alfa * I, -x); //Equivalente a encontrar y = (A - alfa*I)^-1 *x
        lambda = y(p, 1);
        if norm(y, %inf) == 0 then
            disp("A tem autovalor 0.");
            break;
        end
        y = y / y(p, 1);
        err = norm(x - y, %inf);
        x = y;
        if err < epsilon
            disp("Procedimento bem-sucedido.");
            break;
        end
        N = N + 1;
    end
    if N > M then
        N = N - 1;
        disp("Número máximo de iterações excedido.");
    end
    lambda = 1 / lambda + alfa; //O autovalor que encontramos é da matriz (A - alfa*I)^-1. Se lambda é um autovalor de A, então lambda_0 = 1/(lambda - q) é autovalor de (A - alfa*I)^-1. Apenas isolamos lambda.
endfunction

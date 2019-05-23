function [x, C]=Gaussian_Elimination_1(A, b)
C=[A,b];
[n]=size(C,1);
for j=1:(n-1)
 //O pivô está na posição (j,j)
 for i=(j+1):n
//O elemento C(i,j) é o elemento na posição (i,j) of L na decomposição LU de A
 C(i,j)=C(i,j)/C(j,j);
 //Linha i <- Linha i - C(i,j)*Linha j
//Somente os elementos acima da diagonal são computados (aqueles que
//compõem a matrix U)
 C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1);
 end
end

x=zeros(n,1);
// Calcula x, sendo Ux=C(1:n,n+1)
x(n)=C(n,n+1)/C(n,n);
for i=n-1:-1:1
 x(i)=(C(i,n+1)-C(i,i:n)*x(i:n))/C(i,i);
end
C=C(1:n,1:n);
endfunction

function [x, C, trocas]=Gaussian_Elimination_2(A, b)
    C=[A,b];
    [n]=size(C,1);
    trocas = []
    for j=1:(n-1)
        //O pivô está na posição (j,j)
        if C(j,j) == 0
            //Linhas em que o pivô para j não é nulo
            vec = find(C(:,j));
            //k é a linha com a qual vamos trocar
            k = vec(vec > j)(1);
            //Troca de linhas
            D = C(j,:)
            C(j,:) = C(k,:)
            C(k,:) = D
            trocas = [trocas; j, k];
        end
        for i=(j+1):n
            //O elemento C(i,j) é o elemento na posição (i,j) of L na decomposição LU de A
            C(i,j)=C(i,j)/C(j,j);
            //Linha i <- Linha i - C(i,j)*Linha j
            //Somente os elementos acima da diagonal são computados (aqueles que
            //compõem a matrix U)
            C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1);
        end
    end

    x=zeros(n,1);
    // Calcula x, sendo Ux=C(1:n,n+1)
    x(n)=C(n,n+1)/C(n,n);
    for i=n-1:-1:1
        x(i)=(C(i,n+1)-C(i,i:n)*x(i:n))/C(i,i);
    end
    C=C(1:n,1:n);
endfunction

function [x, C, trocas]=Gaussian_Elimination_3(A, b)
    C=[A,b];
    [n]=size(C,1);
    trocas = []
    for j=1:(n-1)
        //O pivô está na posição (j,j)
        if C(j,j) == 0
            //Índice da linha com o maior pivô abaixo de C(j,j))
            [maior, k] = max(abs(C(j+1:n, j)))
            k = k + j;
            //Troca de linhas
            D = C(j,:)
            C(j,:) = C(k,:)
            C(k,:) = D
            trocas = [trocas; j, k];
        end
        for i=(j+1):n
            //O elemento C(i,j) é o elemento na posição (i,j) of L na decomposição LU de A
            C(i,j)=C(i,j)/C(j,j);
            //Linha i <- Linha i - C(i,j)*Linha j
            //Somente os elementos acima da diagonal são computados (aqueles que
            //compõem a matrix U)
            C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1);
        end
    end

    x=zeros(n,1);
    // Calcula x, sendo Ux=C(1:n,n+1)
    x(n)=C(n,n+1)/C(n,n);
    for i=n-1:-1:1
        x(i)=(C(i,n+1)-C(i,i:n)*x(i:n))/C(i,i);
    end
    C=C(1:n,1:n);
endfunction

function [x, C, trocas]=Gaussian_Elimination_4(A, b)
    C=[A,b];
    [n]=size(C,1);
    trocas = []
    for j=1:(n-1)
        //O pivô está na posição (j,j)
        //Índice da linha com o maior pivô abaixo de C(j,j))
        [maior, k] = max(abs(C(j:n, j)))
        k = k + j - 1;
        //Troca de linhas
        if j ~= k
            D = C(j,:)
            C(j,:) = C(k,:)
            C(k,:) = D
            trocas = [trocas; j, k];
        end  
        for i=(j+1):n
            //O elemento C(i,j) é o elemento na posição (i,j) of L na decomposição LU de A
            C(i,j)=C(i,j)/C(j,j);
            //Linha i <- Linha i - C(i,j)*Linha j
            //Somente os elementos acima da diagonal são computados (aqueles que
            //compõem a matrix U)
            C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1);
        end
    end

    x=zeros(n,1);
    // Calcula x, sendo Ux=C(1:n,n+1)
    x(n)=C(n,n+1)/C(n,n);
    for i=n-1:-1:1
        x(i)=(C(i,n+1)-C(i,i:n)*x(i:n))/C(i,i);
    end
    C=C(1:n,1:n);
endfunction

function X = Resolve_com_LU(C, B, trocas)
    [n, m] = size(B);
    //Faz as trocas em B
    for t = 1:size(trocas, 1)
        aux = B(trocas(t, 1),:);
        B(trocas(t, 1),:) = B(trocas(t, 2),:);
        B(trocas(t, 2),:) = aux; 
    end
    //X tem mesma ordem de b
    X = zeros(n, m);
    //Cada sistema:
    for i=1:m
        b = B(:,i);
        //Ax = b
        //(LU)x = b
        //L(Ux) = b
        //Ly = b
        //Calcula y
        y = zeros(n, 1);
        y(1) = b(1); //Seria uma divisão por 1
        for j = 2:n
            y(j) = b(j) - C(j, 1:j - 1) * y(1:j - 1);
        end
        //Ux = y
        //Calcula x
        x=zeros(n, 1);
        // Calcula x, sendo Ux=y
        x(n)=y(n)/C(n,n);
        for k=n-1:-1:1
            x(k)= (y(k) - C(k, k+1:n) * x(k + 1:n)) / C(k, k); 
        end
        X(:, i) = x;
    end
endfunction

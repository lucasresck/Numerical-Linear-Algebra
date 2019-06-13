function A = compression(A, p)
    A = double(A)
    r = rank(A)
    s = max(1, floor(p*r))
    [U, S, V] = svd(A)
    A = U(:, 1:s)*S(1:s, 1:s)*V'(1:s, :)
    A = iconvert(A, 11)    
endfunction

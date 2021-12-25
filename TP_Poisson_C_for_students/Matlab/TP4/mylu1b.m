function [L,U] = mylu1b(A)
    [n, n2] = size(A); 
    for k = [1:n-1]
        A(k+1:n, k) = A(k+1:n, k)/A(k, k);
        A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - A(k+1:n,k) * A(k,k+1:n);
    end 
    % creation L, U pour affichage et verification
    L = A ; 
    U = A ; 
    for i = [1:n]
        L(i,i) = 1 ; 
    end
    for i = [2:n]
        for j = [1:i-1]
            L(j,i) = 0;
            U(i,j) = 0 ;
        end
    end

end

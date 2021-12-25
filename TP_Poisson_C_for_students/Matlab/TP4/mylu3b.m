function [L,U] = mylu3b(A)
    [n, n2] = size(A); 
    for k = [1:n-1]
        for i = [k+1:n]
            A(i, k) = A(i, k)/A(k, k);
        end 
        for i = [k+1:n]
            for j = [k+1:n]
                A(i,j) = A(i,j) - A(i,k) * A(k,j);
            end 
        end 
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

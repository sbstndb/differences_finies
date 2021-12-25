function [L,U] = lutrig(A)
    [n, n2] = size(A); 
    L = eye(n,n);
    U = zeros(n,n)
    U(1,1) = A(1,1) ; 
    L(2,1) = A(2,1) / A(1,1) ;
    U(1,2) = A(1,2);
    for i = [2:n-1]
        U(i,i+1) = A(i,i+1);
        U(i,i) = A(i,i) - A(i-1, i) * L(i, i-1);
        L(i+1, i) = A(i+1, i) / U(i,i);
    end
    U(n,n) =  A(n,n) - A(n-1, n) * L(n, n-1);
end
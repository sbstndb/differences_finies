function [A] = lutrigplace(A)
    [n, n2] = size(A); 
    L = eye(n,n);
    U = zeros(n,n)
    A(1,1) = A(1,1) ; 
    A(1,2) = A(1,2);
    A(2,1) =  A(2, 1) / A(1,1)
    for i = [2:n-1]
        A(i,i+1) = A(i,i+1);
        A(i,i) = A(i,i) - A(i-1, i) * A(i, i-1);
        A(i+1, i) = A(i+1, i) / A(i,i);
    end
    A(n,n) =  A(n,n) - A(n-1, n) * A(n, n-1);
end
function [L,D] = choleskyb(A)
[n1, n2] = size(A); 
L = eye(n1,n1);
D = zeros(n1,n1);
v = zeros(n1,1);
for j = [1:n1]
    for i = [1:j-1]
        v(i) = L(j,i) * D(i,i) ; 
    end
    D(j,j) = A(j,j) - L(j, 1:j-1) * v(1: j-1) ;
    L(j+1:n1, j) = (1/D(j,j))*(A(j+1:n1, j) - L(j+1:n1, 1:j-1)*v(1:j-1)) ; 
end
function [A] = choleskybp(A)
[n1, n2] = size(A); 
v = zeros(n1,1);
A(1+1:n1, 1) = (1/A(1,1))*A(1+1:n1, 1)
for j = [2:n1]
    for i = [1:j-1]
        v(i) = A(j,i) * A(i,i) ; 
    end
    disp(A(j, 1:j-1) * v(1: j-1))
    A(j,j) = A(j,j) - A(j, 1:j-1) * v(1: j-1) ;
    A(j+1:n1, j) = (1/A(j,j))*(A(j+1:n1, j) - A(j+1:n1, 1:j-1)*v(1:j-1)) ; 
end
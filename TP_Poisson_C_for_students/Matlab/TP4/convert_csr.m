function [values, p, s] = convert_csr(A)
[nx, ny] = size(A);
occ = 0 ;
for i = [1:nx]
    for j = [1:ny]
        if A(i,j) ~= 0 
            occ = occ + 1 ; 
        end
    end
end
values = zeros(occ,1) ; 
p = zeros(ny, 1);
s = zeros(nx+1, 1);  % devrait etre resize
pos = 0
size_s = 1 ; 
s(1) = size_s ; 
for i = [1:nx]
    for j = [1:ny]
        if (A(i,j) ~= 0)
            %% valeur non nulle 
            pos = pos + 1 ;
            size_s = size_s + 1 ; 
            values(pos) = A(i,j) ;
            p(pos) = j ; 
        end
    end
    s(i+1) =  size_s ;
end
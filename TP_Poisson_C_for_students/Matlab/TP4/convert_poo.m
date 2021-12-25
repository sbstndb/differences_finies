function [values, px, py] = convert_poo(A)
[nx, ny] = size(A);
occ = 0 ;
for i = [1:nx]
    for j = [1:ny]
        if A(i,j) ~= 0 
            occ = occ + 1 ; 
        end
    end
end
values = zeros(occ, 1);
px = zeros(occ,1);
py = zeros(occ,1);
pos = 0
for i = [1:nx]
    for j = [1:ny]
        if A(i,j) ~= 0 
            pos = pos + 1 ;
            values(pos) = A(i,j) ; 
            px(pos) = i ;
            py(pos) = j ;
        end
    end
end


function [y] = mydspmv(A, CA, RA, x)
    %[n1, n] = size(A);
    n = 8
    y = zeros(n);

    for i = 1 : n
        disp("i")
        disp(i)
        for j = [RA(i):RA(i + 1)-1]
            disp("j")
            disp(j)

            y(i) = y(i) + A(j) * x(CA(j));
        end
    end
end
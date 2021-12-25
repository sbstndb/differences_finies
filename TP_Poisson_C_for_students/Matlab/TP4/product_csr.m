function [v_r] = product_csr(v_a, p_a, s_a, v)
%% _a : for A matric
% v_ : values
% p_ : coordinates
% s_ : coordinate2
[n1, n2] = size(s_a);
n1 = n1 - 1 ; 
v_r = zeros(n1, 0) ; 
index = 1
for i = [1:n1]
        sum = 0  ; 
        for j = [1 : s_a(i+1) - s_a(i)]
            sum = sum + v_a(s_a(i) + j - 1) * v(p_a(index)) ; 
            %disp("ij")
            %disp(v_a(s_a(i) + j - 1))
            %disp(v(p_a(index)))
            index = index + 1 ; 
        end
        v_r(i) = sum ; 
end
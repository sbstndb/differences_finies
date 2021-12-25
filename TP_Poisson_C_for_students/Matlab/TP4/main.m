%%%%%%%%%%% TP4 , exercice 1

size_s = 4 ; 
% generation mtrice symetrique 
A_sym = zeros(size_s,size_s);
for i = [1:size_s]
    for j = [1:size_s]
    A_sym(i,j) = floor(rand()*10);
    end
end

A_sym = A_sym + A_sym' ; 
[L_sym, D_sym] = choleskyb(A_sym) ; 
LT_sym = L_sym' ; 
A_place = choleskybp(A_sym);

% affichage
disp("A_sym")
disp(A_sym)
disp("L")
disp(L_sym)
disp("D")
disp(D_sym)
disp("A en place")
disp(A_place)
disp("L*D*LT")
disp(L_sym*D_sym*LT_sym)

% exemple de matrice symetrique non-definie positive :
A_ndp = [1 0 0; 0 -1 0 ; 0 0 1];
[L_ndp, D_ndp] = choleskyb(A_ndp);



%% mesure de temps 


ech = 5;
power = [1:ech];
time_2 =[];
time_lu_1 = [] ; 
time_lu_3 = [] ; 
time_matlab = []
N_2 = [];
N_lu_1 =[];
N_lu_3 =[];
N_matlab =[];
tic()
min_time = toc();

for i = power
    if (i < 4)
        incm = 10000 ; 
    elseif( i < 7)
        incm = 5000 ; 
    elseif (i < 9)
        incm = 100 ; 
    elseif( i < 11)
        incm = 10;
    else 
        incm = 1 ; 
    end
    n = 2^i;
    disp(i)
    disp(n)
    A = rand(n, n);
    A = A + A' ; %% symmetrisation
    A_p = A ; 
    LLU = A_p ; 
    ULU = A_p;
    if (n < 8000)
        tic()
        for inc = [1: incm]
            %[LLU, ULU] = choleskyb(A);
            choleskyb(A);
        end
        t = toc()/inc;
        time_2 = [time_2, t];
        N_2 = [N_2,n];
    end
    if (n < 4000)
        tic();
        for inc = [1: incm]
            mylu3b(A);
        end
        t = toc()/inc;
        time_lu_3 = [time_lu_3, t];
        N_lu_3 = [N_lu_3,n];
    end
    if (n < 8000)
        tic();
        for inc = [1: incm]
            mylu1b(A);
        end
        t = toc()/inc;
        time_lu_1 = [time_lu_1, t];
        N_lu_1 = [N_lu_1,n];
    end
     if (n < 32000)
        tic();
        for inc = [1: incm]
            lu(A);
        end
        t = toc()/inc;
        time_matlab = [time_matlab, t];
        N_matlab = [N_matlab,n];
    end
    
  
end
loglog(N_lu_3,time_lu_3, '+', color = 'red')
hold on 
loglog(N_lu_1,time_lu_1, '+', color = 'green')
hold on
loglog(N_2,time_2, "+", color = 'blue')
hold on
loglog(N_matlab,time_matlab, "+", color = 'black' )
hold on
loglog(N_lu_3,time_lu_3, color = 'red')
hold on 
loglog(N_lu_1,time_lu_1, color = 'green')
hold on
loglog(N_2,time_2, color = 'blue')
hold on
loglog(N_matlab,time_matlab, color = 'black')
xlabel("size of matrix")
ylabel("time in seconds")
grid()

hold on
legend({"lu3b", "lu1b", "choleskyb", "Matlab"})
legend("Location", "northwest")




%%%%%%%%%%%%%%%%% TP4 , exercice 2
size_s = 4 ; 

% generation 
A_tri = zeros(size_s, size_s);
for i = [1:size_s-1]
    A_tri(i,i) = floor(rand()*10);
    A_tri(i+1, i) = floor(rand()*10) ; 
    A_tri(i, i+1) = floor(rand()*10);
end
A_tri(size_s, size_s) = floor(rand()*10);
[L,U] = lutrig(A_tri) ; 
A_r = L*U ; 
A_p = lutrigplace(A_tri) ;

% diff
diff = norm(A_tri - A_r, "inf");
disp("A")
disp(A_tri)
disp("L")
disp(L)
disp("U")
disp(U)
disp("L*U")
disp(L*U)
disp("Factorisation en place")
disp(A_r)
disp("Difference A - L*U")
disp(diff)



%%%%%%%%%%%%%%% TP4 , exercice 5
A = [15 0 0 22 0 -15 0 0;
0 11 3 0 0 0 2 0;
0 0 0 -6 0 0 0 0;
0 0 0 0 0 0 0 0;
91 0 0 0 0 0 25 7;
0 0 28 0 0 0 0 -2]

[v_poo, px_poo, py_poo] = convert_poo(A);
[v_csr, p_csr, s_csr] = convert_csr(A) ; 
disp("A")
disp(A);
disp("values - POO")
disp(v_poo')
disp("px - POO")
disp(px_poo')
disp("py - POO")
disp(py_poo')
disp("values - CSR")
disp(v_csr')
disp("p - CSR")
disp(p_csr')
disp("s - CSR")
disp(s_csr')


%% produit matrice vecteur
% vecteur 
%v = rand(1,8);
v = [ 1 2 3 4 5 6 7 8];
[v_r] = product_csr(v_csr, p_csr, s_csr, v);
disp("v")
disp(v)
disp("A*v")
disp((A*v')')
disp("product_csr(A_csr, p_csr, s_csr, v)")
disp(v_r)

v = [ 1 1 1 1 1 1 1 1];
[v_r] = product_csr(v_csr, p_csr, s_csr, v);
disp("v")
disp(v)
disp("A*v")
disp((A*v')')
disp("product_csr(A_csr, p_csr, s_csr, v)")
disp(v_r)

v = [ 1 0 0 1 0 0 1 0 ];
[v_r] = product_csr(v_csr, p_csr, s_csr, v);
disp("v")
disp(v)
disp("A*v")
disp((A*v')')
disp("product_csr(A_csr, p_csr, s_csr, v)")
disp(v_r)




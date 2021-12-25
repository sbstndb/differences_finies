%%%%
n = 5 ; % dimension du probleme
T0 = -5.0 ;
T1 = 5.0 ; 


% creation des vecteurs de sauvegarde de la convergence
dimw = 10 ; % nombre de w utilises
dimn = 1 ; % nombre de n utilises 
dima = 20 ; 

iteration_max = 300 ; 
n_iteration = 300;
sauvegarde_j_2 = zeros(dimw, dimn, n_iteration) ; 
sauvegarde_gs_2 = zeros(dimw, dimn, n_iteration) ; 
sauvegarde_r_2 = zeros(dima, dimn, n_iteration) ; 
sauvegarde_j_max = zeros(dimw, dimn, n_iteration) ; 
sauvegarde_gs_max = zeros(dimw, dimn, n_iteration) ; 
sauvegarde_r_max = zeros(dima, dimn, n_iteration) ; 
sauvegarde_w = zeros(dimw,1) ; 
sauvegarde_n = zeros(dimn, 1) ; 
sauvegarde_i = zeros(n_iteration, 1) ;  
for nn = [1:dimn]
    tic();
    disp(" nn  : ") ; 
    disp(nn);
    %n = 2^(nn+2); 
    %n = 2 + nn * 1 + 40*4 ; 
    n = 100 ; 
    sauvegarde_n(nn) = n ;

    A = 2 * eye(n) ;
    A(2:n, 1:n-1) = A(2:n, 1:n-1) -1 * eye(n-1) ; 
    A(1:n-1, 2:n) = A(1:n-1, 2:n) -1 * eye(n-1) ;

    D = 2 * eye(n) ; 

    N = zeros(n) ;
    N(2:n, 1:n-1) = N(2:n, 1:n-1) +1 * eye(n-1) ; 
    N(1:n-1, 2:n) = N(1:n-1, 2:n) +1 * eye(n-1) ; 
    E = tril(N) ; 
    F = triu(N) ; 

    P = D - E ; 
    PinvGS = (D - E)^-1 ; 

    A = sparse(A) ;  
    D = sparse(D) ; 
    N = sparse(N) ;  
    F = sparse(F) ; 
    E = sparse(E) ;     
    P = sparse(P) ; 
    PinvGS = sparse(PinvGS) ; 



    b = zeros(n,1) ;
    b(1) = T0 ; 
    b(n) = T1 ; 
   % b = sparse(b);
    x_direct = A\b ;

    for ww = [1:dimw]
        %w = (ww) / (dimw) ;
        %sauvegarde_w(ww) = 2 *w ;
        w =  (ww) / (dimw) ;
        sauvegarde_w(ww) = w ;

        xj = zeros(n , 1) ; 
        xgs = zeros(n , 1) ;  
        %xj = sparse(xj) ; 
        %xgs = sparse(xgs) ; 

        for iiteration = [1:n_iteration]
            pas_iteration =  iteration_max / n_iteration ;
            sauvegarde_i(iiteration) = iiteration*pas_iteration ;         
            [xj] = jacobim(0.5, N, b, xj, w, pas_iteration) ;   % methode de jacobi avec sur relaxation 
            [xgs] = gauss_seidelm(PinvGS, F , b, xgs, w, pas_iteration); ; % methode de gauss-seidel avec sur relaxation

            % calcul de l'erreur 
            errj = norm(xj - x_direct,2)/norm(x_direct,2) ; 
            errgs = norm(xgs - x_direct,2)/norm(x_direct,2) ; 
            sauvegarde_j_2(ww, nn, iiteration) = errj ; 
            sauvegarde_gs_2(ww, nn, iiteration) = errgs ;          
            errj = norm(xj - x_direct,inf)/norm(x_direct,inf) ; 
            errgs = norm(xgs - x_direct,inf)/norm(x_direct,inf) ; 
            sauvegarde_j_max(ww, nn, iiteration) = errj ; 
            sauvegarde_gs_max(ww, nn, iiteration) = errgs ;   
        end
    end
 for aa = [1:dima]
     disp("aa")
     disp(aa)
        %w = (ww) / (dimw) ;
        %sauvegarde_w(ww) = 2 *w ;
        a =   0.50005 * (aa) / (dima) ;
        sauvegarde_a(aa) = a ;

        xj = zeros(n , 1) ; 
        xgs = zeros(n , 1) ;  
        xa = zeros(n, 1) ; 
        %xj = sparse(xj) ; 
        %xgs = sparse(xgs) ; 

        for iiteration = [1:n_iteration]
            pas_iteration =  iteration_max / n_iteration ;
            sauvegarde_i(iiteration) = iiteration*pas_iteration ;         
            [xa] = richarson(sparse(eye(n,n)), A, b, xa, a, pas_iteration) ;
            % calcul de l'erreur 
            erra = norm(xa - x_direct,2)/norm(x_direct,2) ; 
            sauvegarde_r_2(aa, nn, iiteration) = erra ;     
            erra = norm(xa - x_direct,inf)/norm(x_direct,inf) ; 
            sauvegarde_r_max(aa, nn, iiteration) = erra ; 
        end
    end

    toc();
end

%% test avec les formulations matricielles 
%xj = zeros(n,1);
%xgs = zeros(n,1);
%[xj] = jacobim(0.5, N, b, xj, 1, 10 * pas_iteration);
%[xgs] = gauss_seidelm(PinvGS, F , b, xj, 0.5, 10 * pas_iteration);



%% plot 
% plot de la convergence des methodes en fonction de la taille du probleme
% - avec w = 1.0

% jacobi, w = 1, norme 2
% tmp = zeros(n_iteration, dimn) ; 
% for i = [1:n_iteration] 
%     tmp(i,:) = sauvegarde_j_2(dimw,:,i) ; 
% end
% 
% for i  = [1:dimn]
%     tmp2 = tmp(:,i) ; 
%     p3 = plot3(sauvegarde_i, sauvegarde_n(i)*ones(n_iteration,1), tmp2)
%     hold on ;
% end
% set(gca, 'zscale', 'log')
% set(gca, 'yscale', 'log')
% view(40,20)
% xlabel("iteration")
% ylabel("taille n de matrice")
% zlabel("erreur relative norme 2")
% grid on
% savefig("j_i_n_2.fig")
% close()

% gs, w = 1, norme 2
% tmp = zeros(n_iteration, dimn) ; 
% for i = [1:n_iteration] 
%     tmp(i,:) = sauvegarde_gs_2(dimw,:,i) ; 
% end
% 
% for i  = [1:dimn]
%     tmp2 = tmp(:,i) ; 
%     p3 = plot3(sauvegarde_i, sauvegarde_n(i)*ones(n_iteration,1), tmp2)
%     hold on ;
% end
% set(gca, 'zscale', 'log')
% set(gca, 'yscale', 'log')
% view(40,20)
% xlabel("iteration")
% ylabel("taille n de matrice")
% zlabel("erreur relative norme 2")
% grid on
% savefig("gs_i_n_2.fig")
% close()

% jacobi, n = max, w variable, norme 2
% tmp = zeros(n_iteration, dimw) ; 
% for i = [1:n_iteration] 
%    tmp(i,:) = sauvegarde_j_2(:,dimn,i)' ; 
% end
% 
% for i  = [1:dimw]
%    
%    tmp2 = tmp(:,i) ; 
%    p3 = plot3(sauvegarde_i, sauvegarde_w(i)*ones(n_iteration,1), tmp2)
%    hold on ;
% end
% set(gca, 'zscale', 'log')
% set(gca, 'yscale', 'log')
% view(40,20)
% xlabel("iteration")
% ylabel("valeur de w")
% zlabel("erreur relative norme 2")
% grid on
% savefig("j_i_w_2.fig")
% close()
% 
% % jacobi, n = max, w variable, norme 2
% tmp = zeros(n_iteration, dimw) ; 
% for i = [1:n_iteration] 
%    tmp(i,:) = sauvegarde_gs_2(:,dimn,i)' ; 
% end
% 
% for i  = [1:dimw]
%    tmp2 = tmp(:,i) ; 
%    p3 = plot3(sauvegarde_i, sauvegarde_w(i)*ones(n_iteration,1), tmp2)
%    hold on ;
% end
% set(gca, 'zscale', 'log')
% set(gca, 'yscale', 'log')
% view(40,20)
% xlabel("iteration")
% ylabel("valeur de w")
% zlabel("erreur relative norme 2")
% grid on
% savefig("gs_i_w_2.fig")
% close()


%%% norme inf 
% tmp = zeros(n_iteration, dimn) ; 
% for i = [1:n_iteration] 
%     tmp(i,:) = sauvegarde_j_max(dimw,:,i) ; 
% end
% 
% for i  = [1:dimn]
%     tmp2 = tmp(:,i) ; 
%     p3 = plot3(sauvegarde_i, sauvegarde_n(i)*ones(n_iteration,1), tmp2)
%     hold on ;
% end
% set(gca, 'zscale', 'log')
% set(gca, 'yscale', 'log')
% view(40,20)
% xlabel("iteration")
% ylabel("taille n de matrice")
% zlabel("erreur relative norme inf")
% grid on
% savefig("j_i_n_inf.fig")
% close()

% gs, w = 1, norme 2
% tmp = zeros(n_iteration, dimn) ; 
% for i = [1:n_iteration] 
%     tmp(i,:) = sauvegarde_gs_max(dimw,:,i) ; 
% end
% 
% for i  = [1:dimn]
%     tmp2 = tmp(:,i) ; 
%     p3 = plot3(sauvegarde_i, sauvegarde_n(i)*ones(n_iteration,1), tmp2)
%     hold on ;
% end
% set(gca, 'zscale', 'log')
% set(gca, 'yscale', 'log')
% view(40,20)
% xlabel("iteration")
% ylabel("taille n de matrice")
% zlabel("erreur relative norme inf")
% grid on
% savefig("gs_i_n_inf.fig")
% close()

% jacobi, n = max, w variable, norme 2
% tmp = zeros(n_iteration, dimw) ; 
% for i = [1:n_iteration] 
%     tmp(i,:) = sauvegarde_j_max(:,dimn,i)' ; 
% end
% 
% for i  = [1:dimw]
%     
%     tmp2 = tmp(:,i) ; 
%     p3 = plot3(sauvegarde_i, sauvegarde_w(i)*ones(n_iteration,1), tmp2)
%     hold on ;
% end
% set(gca, 'zscale', 'log')
% set(gca, 'yscale', 'log')
% view(40,20)
% xlabel("iteration")
% ylabel("valeur de w")
% zlabel("erreur relative norme inf")
% grid on
% savefig("j_i_w_inf.fig")
% close()
% 
% %jacobi, n = max, w variable, norme 2
% tmp = zeros(n_iteration, dimw) ; 
% for i = [1:n_iteration] 
%     tmp(i,:) = sauvegarde_gs_max(:,dimn,i)' ; 
% end
% 
% for i  = [1:dimw]
%     tmp2 = tmp(:,i) ; 
%     p3 = plot3(sauvegarde_i, sauvegarde_w(i)*ones(n_iteration,1), tmp2)
%     hold on ;
% end
% set(gca, 'zscale', 'log')
% set(gca, 'yscale', 'log')
% view(40,20)
% xlabel("iteration")
% ylabel("valeur de w")
% zlabel("erreur relative norme inf")
% grid on
% savefig("gs_i_w_inf.fig")
% close()



%% post traitement de richardson
maxindex = 0 ; 
for i = [1:dima]
       if (sauvegarde_r_2(i,1, n_iteration) < 0.999999999)
           maxindex = i ; 
       end
end

maxalpha = sauvegarde_a(maxindex) ; 
disp("a maximun avant divergence : ")
disp(maxalpha)

% for n = 1000 : alpha lim between  0.5170 nd 0.5180
% for n = 10 : alpha lim between  0.5180 nd 0.5190


%richardson, n = max, a variable, norme 2
tmp = zeros(n_iteration, dima) ; 
for i = [1:n_iteration] 
    tmp(i,:) = sauvegarde_r_2(:,dimn,i)' ; 
end

for i  = [1:dima]
    
    tmp2 = tmp(:,i) ; 
    p3 = plot3(sauvegarde_i, sauvegarde_a(i)*ones(n_iteration,1), tmp2) ;
    hold on ;
end
set(gca, 'zscale', 'log')
%set(gca, 'yscale', 'log')
view(70,20)
xlabel("iteration")
ylabel("valeur de a")
zlabel("erreur relative norme 2")
grid on
savefig("j_i_a_2.fig")
close()


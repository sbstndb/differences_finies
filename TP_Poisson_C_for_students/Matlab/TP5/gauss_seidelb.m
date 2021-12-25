function [x] = gauss_seidelb(A, b, xinit, w, iteration)
%% A est la matrice du probleme
%% b est le second membre
%% xinit est le vecteur solution a l'etat initial
%% iteration est le nombre d'iterations du schema iteratif
%% sortie : x est le vecteur solution
xold = xinit ;
[n,n1] = size(b) ;
for it = [1 : iteration]
    xnew = zeros(n,1) ; 
    for j = [1:n]
        xnew(j) = b(j) ;
        for k = [1:j-1]
            xnew(j) = xnew(j) - A(j,k) * xnew(k) ;
        end
        for k = [j+1:n]
            xnew(j) = xnew(j) - A(j,k) * xold(k) ;
        end   
        xnew(j) = w * xnew(j) / A(j,j) + (1 - w) * xold(j) ; 
    end
    xold = xnew ; 
end
x = xnew ; 
end

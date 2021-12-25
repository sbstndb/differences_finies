function [x] = jacobim(Dinv, N, b, x, w, iteration)
%% Dinv est la valeur diagonale inverse de A
%% N est N = A - D(A), soit A privé de sa diagonale
%% b est le second membre
%% xinit est le vecteur solution a l'etat initial
%% iteration est le nombre d'iterations du schema iteratif
%% sortie : x est le vecteur solution
[n,n1] = size(b) ;
for it = [1 : iteration]
    x = w * Dinv * ( N * x + b) + (1 - w) * x ;  
end
end

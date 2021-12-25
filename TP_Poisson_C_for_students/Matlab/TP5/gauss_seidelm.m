function [x] = gauss_seidelm(Pinv, N, b , x, w, iteration)
%% A est la matrice du probleme
%% b est le second membre
%% xinit est le vecteur solution a l'etat initial
%% iteration est le nombre d'iterations du schema iteratif
%% sortie : x est le vecteur solution
[n,n1] = size(b) ;
for it = [1 : iteration]
    x = w * Pinv * ( N * x + b) + (1 - w) * x ; 
end
end

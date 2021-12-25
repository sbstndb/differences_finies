function [x] = richardson(Pinv, A, b , x, alpha, iteration)
%% Pinv est la matrice de preconditionnement inversé 
%% A est la matrice du probleme
%% b est le second membre
%% x est le vecteur solution a l'etat initial
%% iteration est le nombre d'iterations du schema iteratif
%% sortie : x est le vecteur solution
[n,n1] = size(b) ;
for it = [1 : iteration]
    x = x + alpha * Pinv * (b - A * x) ; 
end
end

function G = loadGauntMatrices(lmax, d)
% 
% function to load pre-computed Gaunt matrices of order d.
% 
% lmax is the desired max degree of the Laplace expasion.
% d is the dimension of the Gaunt matrix
% 

global G3 G4 G5

if lmax > 8
    error('AGaun matrices for values of lmax higher thab 8 have not been precomputed!!');
end
if d < 3 || d > 5
    error('minimum value for d is 3, maximum values is 5');
end


switch d
    case 3
        load(G3, 'G');
    case 4
        load(G4, 'G');
    case 5
        load(G5, 'G');
end


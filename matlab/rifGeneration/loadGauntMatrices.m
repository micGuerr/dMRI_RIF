function G = loadGauntMatrices(lmax, d)
% 
% function to load pre-computed Gaunt matrices of order d.
% 
% lmax is the desired max degree of the Laplace expasion.
% d is the dimension of the Gaunt matrix
% 

loadFunPath = which('loadGauntMatrices');
setPath = fileparts(loadFunPath);

if lmax > 8
    error('AGaun matrices for values of lmax higher thab 8 have not been precomputed!!');
end

switch d
    case 3
        gPath = fullfile(setPath, 'Gd3'); 
    case 4
        gPath = fullfile(setPath, 'Gd4'); 
    case 5
        gPath = fullfile(setPath, 'Gd5'); 
end

load(gPath, 'G');
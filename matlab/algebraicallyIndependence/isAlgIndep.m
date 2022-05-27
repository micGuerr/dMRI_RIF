function a = isAlgIndep(X,P)
% 
% Check if a set of ploynomials is alebraically independent
% 
% Usage:
%   a = isAlgIndep(X,P)
% 
% Input:
%   X
% 
% 
% 
% 
% 

nX = length(X);
nP = length(P);

J = computeJacobian(X,P);

% generate 100 different random values of X and check if J has full rank
N = 100;
x = rand(N, nX);
jrank = nan(N,1);

for i = 1:N
    jval = double(subs(J, X(1:end), x(i,1:end)));
    jrank(i) = rank(jval);
end
    
a = sum(jrank>=nP) == N;



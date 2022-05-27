function J = computeJacobian(X,P)
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 


nX = length(X);
nP = length(P);

J = sym('J', [nP nX]);

for pp = 1:nP
    for xx = 1:nX
        J(pp,xx) = diff(P(pp), X(xx));
    end
end
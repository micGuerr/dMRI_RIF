function U = Im2Real_SH(l)
% 
% Given a specific value of l it computes the (2l+1)x(2l+1) unitary matrix
% to go from an imaginary set of SH to a real one.
% The computation is based on https://doi.org/10.1016/S0166-1280(96)90531-X
% Note That the signs in the matrix are different to match Desconteaux
% definition of real SH.
% 
% Usage:
%   function U = Im2Real_SH(l)
%
% Input:
%   l   is a natural number 
%
% Output:
%   U   (2l+1)x(2l+1) unitary matrix. If m and mu are vectors containing
%       values from -l to l, m = mu = [-l -l+1 ... l], r and c are the 
%       indices for the rows and columns respct, then the values of U are:
%       if r ~= c
%           U(r,c) =              0
%       else
%           U(r,c) =            ( 1/sqrt(2))    if mu(r) > 0 and m(c) > 0
%           U(r,c) = (-1)^m(r) *( 1/sqrt(2))    if mu(r) > 0 and m(c) < 0
%           U(r,c) =            ( i/sqrt(2))    if mu(r) < 0 and m(c) > 0
%           U(r,c) = (-1)^m(r) *(-i/sqrt(2))    if mu(r) < 0 and m(c) < 0
%           U(r,c) =              1             if mu(r) = 0 and m(c) = 0
%
% Author:
%   Michele Guerreri [m.guerreri@ucl.ac.uk]
%

% Define the m values based on l
m = -l:1:l;
mu = m;

% Define U
U = zeros((2*l+1), (2*l+1));

% Fill the U matrix
for rr = 1:(2*l+1) % loop over all possible mu values
    for cc = [rr (2*l+2)-rr]
        % All off-diagonal values are 0
        if mu(rr) > 0 && m(cc) > 0
            U(rr,cc) = 1/sqrt(2);
        elseif mu(rr) > 0 && m(cc) < 0
            U(rr,cc) = ( (-1)^(m(cc)) )*( 1/sqrt(2) );
        elseif mu(rr) < 0 && m(cc) > 0
            U(rr,cc) = sqrt(-1)/sqrt(2);
        elseif mu(rr) < 0 && m(cc) < 0
            U(rr,cc) = ( -(-1)^(m(cc)) )*( sqrt(-1)/sqrt(2) );
        else
            U(rr,cc) = 1;
        end
    end
end
    





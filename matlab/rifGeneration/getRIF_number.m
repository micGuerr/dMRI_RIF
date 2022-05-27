function n_rif = getRIF_number(L,d)
% 
% Compute the number of rotation invariants feature given a set of values L
% and d. he computation is not analytical but it is based on the table in
% Zucchelli's paper 2018 Neuroimage.
% 
% Usage:
%   n_rif = getRIF_number(L,d)
% 
% Input:
%   L   Nx1 vector of max l degree used for SH rxpansion.
%   d   Nx1 vector of desired degree of RIF.
% 
% Output:
%   n_rif NXi vector with nuber of RIF for each entry
% 
% Author:
%   Michele Guerrei (m.guereri@ucl.ac.uk)
% 
% 

% check that the vectors have same length
N = length(L);
if length(d) ~= N
    error('Dimension mismatch!')
end

% define output
n_rif = nan(N,1);

% loop over the values of L and d
for i = 1:N
    switch d(i)
        case 1
            n_rif(i) = 1;
        case 2
            n_rif(i) = L(i)/2+1;
        case 3
            switch L(i)
                case 2
                    n_rif(i) = 3;
                case 4
                    n_rif(i) = 7;
                case 6
                    n_rif(i) = 13;
                case 8
                    n_rif(i) = 35;
                otherwise
                    error('not implemented for L higher than 8!!')
            end
        case 4
            switch L(i)
                case 2
                    n_rif(i) = 4;
                case 4
                    n_rif(i) = 12;
                case 6
                    n_rif(i) = 28;
                case 8
                    n_rif(i) = 70;
                otherwise
                    error('not implemented for L higher than 8!!')
            end
        case 5
            switch L(i)
                case 2
                    n_rif(i) = 5;
                case 4
                    n_rif(i) = 18;
                case 6
                    n_rif(i) = 49;
                case 8
                    n_rif(i) = 126;
                otherwise
                    error('not implemented for L higher than 8!!')
            end
        otherwise
            error('d is too high!!')
    end
    
end


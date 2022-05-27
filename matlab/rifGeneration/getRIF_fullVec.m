function rif = getRIF_fullVec(SHcoeff, lmax, d)
% 
% Compute a set of rotation invariant features (RIF) of a real even spherical 
% function expanded in a Laplace series using the real set of spherical
% harmonics (SH).
% It uses the approach described in:
% https://doi.org/10.1016/j.media.2019.101597
% 
% Usage:
%   rif = getRIF(SHcoeff, lmax, d)
% 
% Input:
%   rif         output set of RIF.
%
% Output:
%   SHcoeff     coefficients of the Laplace series. The j-th element of the
%               array is expected to be the coefficient with index 
%               j = l^2 + l + m.
%   lmax        maximum degree of the expansion to be used in the RIF 
%               estimation. It should be a positive even number. 
%   d           power of the RIF. Integer between 1 and 3.
% 
% 
% Author:
%   Michele Guerreri [m.guerreri@ucl.ac.uk]
%
% 

%% Check inputs are sensible

% make sure lmax is even
if rem(lmax,2) ~= 0
    error('The maximum l to be used for the RIF evaluation is expected to be even.');
end

% Number of coefficients must be compatible with the selected lmax
N = (lmax+1)*(lmax+2)/2;
if numel(SHcoeff) < N
    error('Based on the input lmax, the number of SH coefficients should be at least equal to %d.', N);
end


% make sure d is integer and max eq to 3
if d ~= round(d)
    error('The power d should be an integer.');
elseif d > 5
    error('Computation of RIFs with a power d up to 5 have been implemented so far.');
end
%% Compute the RIFs

% The number of possible l values is
nL =  lmax/2+1;
% the max number of RIFs is given by n multichoose k due to the ineherent symmetries:
n_rif = factorial(nL+d-1)/( factorial(d) * factorial(nL-1) );
% the multisubsets are:
mss = nmultichoosek(0:2:lmax, d);
% We allow for simultaneous computation of multiple voxels
nVox = size(SHcoeff,1);


switch d
    case 1
        % This corresponds to direction-averaging
        rif = sqrt(4*pi)*SHcoeff(:,1);
    case 2
        % define the output
        rif = zeros(nVox, n_rif);
        % Kroneker's delta
        kdelta = @(x,y) double(eq(x,y));
        % loop over all the combinations of degrees l
        for ii = 1:n_rif
            tmp_rif = 0;
            l1 = mss(ii,1);
            l2 = mss(ii,2);
            % sum over the different m values
            for m1 = -l1:l1
                j1 = (l1^2)/2 + l1/2 + m1 + 1;
                for m2 = -l2:l2
                    j2 = (l2^2)/2 + l2/2 + m2 + 1;
                    tmp_rif = tmp_rif + ...
                        ( SHcoeff(:,j1).*SHcoeff(:,j2).*kdelta(l1,l2)*kdelta(m1,m2) );
                end
            end
            % assigne the result to the output
            rif(:, ii) = tmp_rif;
        end
    case 3
        rif = zeros(nVox, n_rif);
        % Use pre computed Gaunt matrix
        G = loadGauntMatrices(lmax, d);
        %G = real_gaunt_mtx(lmax,lmax, lmax);
        % loop over all the combinations of degrees l
        for ii = 1:n_rif
            tmp_rif = 0;
            l1 = mss(ii,1);
            l2 = mss(ii,2);
            l3 = mss(ii,3);
            % sum over the different m values
            for m1 = -l1:l1
                j1 = (l1^2)/2 + l1/2 + m1 + 1;
                k1 = (l1^2) + l1 + m1 + 1;
                for m2 = -l2:l2
                    j2 = (l2^2)/2 + l2/2 + m2 + 1;
                    k2 = (l2^2) + l2 + m2 + 1;
                    for m3 = -l3:l3
                        j3 = (l3^2)/2 + l3/2 + m3 + 1;
                        k3 = (l3^2) + l3 + m3 + 1;
                        tmp_rif = tmp_rif + ...
                            ( SHcoeff(:,j1).*SHcoeff(:,j2).*SHcoeff(:,j3).*G(k1, k2, k3) );
                    end
                end
            end
            % assigne the result to the output
            rif(:,ii) = tmp_rif;
        end
    case 4
        rif = zeros(nVox, n_rif);
        % Use pre computed Gaunt matrix
        % G = genReal_gaunt_mtx(lmax, d);
        G = loadGauntMatrices(lmax, d);
        % loop over all the combinations of degrees l
        for ii = 1:n_rif
            tmp_rif = 0;
            l1 = mss(ii,1);
            l2 = mss(ii,2);
            l3 = mss(ii,3);
            l4 = mss(ii,4);
            % sum over the different m values
            for m1 = -l1:l1
                j1 = (l1^2)/2 + l1/2 + m1 + 1;
                k1 = (l1^2) + l1 + m1 + 1;
                for m2 = -l2:l2
                    j2 = (l2^2)/2 + l2/2 + m2 + 1;
                    k2 = (l2^2) + l2 + m2 + 1;
                    for m3 = -l3:l3
                        j3 = (l3^2)/2 + l3/2 + m3 + 1;
                        k3 = (l3^2) + l3 + m3 + 1;
                        for m4 = -l4:l4
                            j4 = (l4^2)/2 + l4/2 + m4 + 1;
                            k4 = (l4^2) + l4 + m4 + 1;
                            tmp_rif = tmp_rif + ...
                                ( SHcoeff(:,j1).*SHcoeff(:,j2).*SHcoeff(:,j3).*SHcoeff(:,j4).*G(k1, k2, k3, k4) );
                        end
                    end
                end
            end
            % assigne the result to the output
            rif(:,ii) = tmp_rif;
        end
    case 5
        rif = zeros(nVox, n_rif);
        % Use pre computed Gaunt matrix
        % G = genReal_gaunt_mtx(lmax, d);
        G = loadGauntMatrices(lmax, d);
        % loop over all the combinations of degrees l
        for ii = 1:n_rif
            tmp_rif = 0;
            l1 = mss(ii,1);
            l2 = mss(ii,2);
            l3 = mss(ii,3);
            l4 = mss(ii,4);
            l5 = mss(ii,5);
            % sum over the different m values
            for m1 = -l1:l1
                j1 = (l1^2)/2 + l1/2 + m1 + 1;
                k1 = (l1^2) + l1 + m1 + 1;
                for m2 = -l2:l2
                    j2 = (l2^2)/2 + l2/2 + m2 + 1;
                    k2 = (l2^2) + l2 + m2 + 1;
                    for m3 = -l3:l3
                        j3 = (l3^2)/2 + l3/2 + m3 + 1;
                        k3 = (l3^2) + l3 + m3 + 1;
                        for m4 = -l4:l4
                            j4 = (l4^2)/2 + l4/2 + m4 + 1;
                            k4 = (l4^2) + l4 + m4 + 1;
                            for m5 = -l5:l5
                                j5 = (l5^2)/2 + l5/2 + m5 + 1;
                                k5 = (l5^2) + l5 + m5 + 1;
                                
                                tmp_rif = tmp_rif + ...
                                    ( SHcoeff(:,j1).*SHcoeff(:,j2).*SHcoeff(:,j3).*SHcoeff(:,j4).*SHcoeff(:,j5).* ...
                                    G(k1, k2, k3, k4, k5) );
                            end
                        end
                    end
                end
            end
            % assigne the result to the output
            rif(:,ii) = tmp_rif;
        end
end

function combs = nmultichoosek(values, k)
% Return number of multisubsets or actual multisubsets.
if numel(values)==1 
    n = values;
    combs = nchoosek(n+k-1,k);
else
    n = numel(values);
    combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
    combs = reshape(values(combs),[],k);
end
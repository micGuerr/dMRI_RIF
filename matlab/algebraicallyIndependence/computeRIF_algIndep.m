function [algIndpRif,rifLabels] = computeRIF_algIndep(Lmax,dmax,output_name)
%
% Compute the sequence of independent RIF for each combination of L and d
% given the specified Lmax and dmax.
%
% Usage:
%   computeRIF_algIndep(Lmax,dmax,output_name)
%
% Input:
%   Lmax
%   dmax
%   output_name     specify output name
%
% Output:
%   for each pair of L and d the function outputs an array of logic w/
%   length the maximum length available for that pair. Ones correspond to the
%   independent RIF. The arrays are stored and saved in a cell array of
%   dimensions [ Lmax/2 d ] which is saved based on the "output_name"
%   variable.
%
% Author:
%   Michele Guerreri [m.guerreri@ucl.ac.uk]
%
%

% define the L and d
L = 2:2:Lmax;
D = 1:dmax;

% first define the output cell array
nL = length(L);
nD = length(D);
algIndpRif = cell(nL,nD);
rifLabels = cell(nL,nD);

% Now loop first over the L
cl = 0; % cont the ls
for l = L
    cl = cl+1;
    % for this l compute how many RIF we are going to check the independ
    % The number of possible l values is
    nL =  l/2+1;
    % the max number of RIFs is given by n multichoose k due to the ineherent symmetries:
    tmp_nrif = factorial(nL+D-1)./( factorial(D) * factorial(nL-1) );
    % first n rif cannot be but 1
    tmp_nrif(1) = 1;
    nrif = cumsum(tmp_nrif);
    
    % define the maximum number of SH coefficients for this L
    nSHcoeff = (l+1)*(l+2)/2;
    % generate the coeffiecients which represent the polynomial variables
    X = sym('x', [1 nSHcoeff]);
    
    % create the array of logical indicating independence.
    indp_idx = true(nrif(end),1);
    % create array of cells to save the index label
    idx_labels = cell(nrif(end),1);
    % define an array of functions for this specific d given by ncomb.
    P = sym('p', [1 nrif(end)]);
    % loop over the d
    for d = D
        %% d = 1
        if d == 1
            % now compute the polynomial for d = 1.
            P(1) = sqrt(4*sym(pi))*X(1);
            % assign the first label
            idx_labels{1} = 0;
            
            % store index and labels into the output cell array
            algIndpRif{cl,1} = indp_idx(1);
            rifLabels{cl,1} = idx_labels(1);
        else
            % first check unique combinations of indices
            comb = nmultichoosek(0:2:l, d);
            ncomb = size(comb,1);
            idx_labels(nrif(d-1)+1:nrif(d)) = num2cell(comb,2);
            
            %% d = 2
            if d == 2
                % Kroneker's delta
                kdelta = @(x,y) double(eq(x,y));
                % loop over all the combinations of degrees l
                for ii = 1:ncomb
                    % the set of coefficients I am going to sum
                    M = sym('f', [1 prod(comb(ii,:)*2+1)]);
                    l1 = comb(ii,1);
                    l2 = comb(ii,2);
                    % sum over the different m values
                    cnt = 0;
                    for m1 = -l1:l1
                        j1 = (l1^2)/2 + l1/2 + m1 + 1;
                        for m2 = -l2:l2
                            cnt = cnt+1;
                            j2 = (l2^2)/2 + l2/2 + m2 + 1;
                            M(cnt) = ( X(j1)*X(j2)*kdelta(l1,l2)*kdelta(m1,m2) );
                        end
                    end
                    % assigne the result to the output
                    P(nrif(d-1)+ii) = sum(M);
                    tmpP = P(1:nrif(d-1)+ii);
                    tmpP(indp_idx(1:nrif(d-1)+ii) == 0 ) = [];
                    indp_idx(nrif(d-1)+ii) = isAlgIndep(X,tmpP);
                    %indp_idx(nrif(d-1)+ii) = isAlgIndep(X,P(1:nrif(d-1)+ii));
                end
            end
            %% d = 3
            if d == 3
                % Use pre computed Gaunt matrix
                G = loadGauntMatrices(l, d);
                % loop over all the combinations of degrees l
                for ii = 1:ncomb
                    % the set of coefficients I am going to sum
                    M = sym('f', [1 prod(comb(ii,:)*2+1)]);
                    l1 = comb(ii,1);
                    l2 = comb(ii,2);
                    l3 = comb(ii,3);
                    % sum over the different m values
                    cnt = 0;
                    % sum over the different m values
                    for m1 = -l1:l1
                        j1 = (l1^2)/2 + l1/2 + m1 + 1;
                        k1 = (l1^2) + l1 + m1 + 1;
                        for m2 = -l2:l2
                            j2 = (l2^2)/2 + l2/2 + m2 + 1;
                            k2 = (l2^2) + l2 + m2 + 1;
                            for m3 = -l3:l3
                                cnt = cnt+1;
                                j3 = (l3^2)/2 + l3/2 + m3 + 1;
                                k3 = (l3^2) + l3 + m3 + 1;
                                M(cnt) = ( X(j1)*X(j2)*X(j3)*G(k1, k2, k3) );
                            end
                        end
                    end
                    % assigne the result to the output
                    P(nrif(d-1)+ii) = sum(M);
                    tmpP = P(1:nrif(d-1)+ii);
                    tmpP(indp_idx(1:nrif(d-1)+ii) == 0 ) = [];
                    indp_idx(nrif(d-1)+ii) = isAlgIndep(X,tmpP);
                    %indp_idx(nrif(d-1)+ii) = isAlgIndep(X,P(1:nrif(d-1)+ii));
                end
            end
            %% d = 4
            if d == 4
                % Use pre computed Gaunt matrix
                G = loadGauntMatrices(l, d);
                % loop over all the combinations of degrees l
                for ii = 1:ncomb
                    % the set of coefficients I am going to sum
                    M = sym('f', [1 prod(comb(ii,:)*2+1)]);
                    l1 = comb(ii,1);
                    l2 = comb(ii,2);
                    l3 = comb(ii,3);
                    l4 = comb(ii,4);
                    % sum over the different m values
                    cnt = 0;
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
                                    cnt = cnt+1;
                                    j4 = (l4^2)/2 + l4/2 + m4 + 1;
                                    k4 = (l4^2) + l4 + m4 + 1;
                                    M(cnt) = ( X(j1)*X(j2)*X(j3)*X(j4)*G(k1, k2, k3, k4) );
                                end
                            end
                        end
                    end
                    % assigne the result to the output
                    P(nrif(d-1)+ii) = sum(M);
                    tmpP = P(1:nrif(d-1)+ii);
                    tmpP(indp_idx(1:nrif(d-1)+ii) == 0 ) = [];
                    indp_idx(nrif(d-1)+ii) = isAlgIndep(X,tmpP);
                    %indp_idx(nrif(d-1)+ii) = isAlgIndep(X,P(1:nrif(d-1)+ii));
                end
            end
            %% d = 5
            if d == 5
                % Use pre computed Gaunt matrix
                G = loadGauntMatrices(l, d);
                % loop over all the combinations of degrees l
                for ii = 1:ncomb
                    % the set of coefficients I am going to sum
                    M = sym('f', [1 prod(comb(ii,:)*2+1)]);
                    l1 = comb(ii,1);
                    l2 = comb(ii,2);
                    l3 = comb(ii,3);
                    l4 = comb(ii,4);
                    l5 = comb(ii,5);
                    % sum over the different m values
                    cnt = 0;
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
                                        cnt = cnt+1;
                                        j5 = (l5^2)/2 + l5/2 + m5 + 1;
                                        k5 = (l5^2) + l5 + m5 + 1;
                                        M(cnt) = ( X(j1)*X(j2)*X(j3)*X(j4)*X(j5)*G(k1, k2, k3, k4, k5) );
                                    end
                                end
                            end
                        end
                    end
                    % assigne the result to the output
                    P(nrif(d-1)+ii) = sum(M);
                    tmpP = P(1:nrif(d-1)+ii);
                    tmpP(indp_idx(1:nrif(d-1)+ii) == 0 ) = [];
                    indp_idx(nrif(d-1)+ii) = isAlgIndep(X,tmpP);
                    %indp_idx(nrif(d-1)+ii) = isAlgIndep(X,P(1:nrif(d-1)+ii));
                end
            end
            algIndpRif{cl,d} = indp_idx(1:nrif(d));
            rifLabels{cl,d} = idx_labels(1:nrif(d));
        end
        save(output_name, 'algIndpRif', 'rifLabels')
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
function A = genReal_gaunt_mtx(N, d)
%GENREAL_GAUNT_MTX Construct a matrix of Gaunt coefficients
%
% GENREAL_GAUNT_MTX constructs the d dimensional matrix of Generalized Real
% Gaunt (R-Gaunt) coefficients which represent the instegral of d
% spherical harmonics such as
%
% The Generalized R-Gaunt coefficients are calculated as the product of
% R-Gaunt coefficients as descirbed by equation (40) in
% https://doi.org/10.1016/S0166-1280(96)90531-X
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Michele Guerreri (m.guerreri@ucl.ac.uk)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure d is integer and max eq to 5
if d ~= round(d)
    error('input d should be an integer.');
elseif d > 5
    error('Computation of generalized R-Gaunt coeff with d above 5 have not been implemented yet.');
elseif d < 3
    error('It doesn''t make sense to use d lower than 3...');
end

% Get the R-Gaunt matrix
%G = real_gaunt_mtx(N,N,N);
G = loadGauntMatrices(N, 3);

% Compute generalized GAUNT coeff:

switch d
    case 3
        A = G;
    case 4
        A = zeros((N+1)^2, (N+1)^2, (N+1)^2, (N+1)^2);
        for n = 0:N
            for mu = -n:n
                % Update indices in A and U mat
                q = n*(n+1)+mu;
                
                for n1 = 0:N
                    for mu1 = -n1:n1
                        % Update indices in A and U1 mat
                        q1 = n1*(n1+1)+mu1;
                        
                        for n2 = 0:N
                            for mu2 = -n2:n2
                                % Update indices in A and U2 mat
                                q2 = n2*(n2+1)+mu2;
                                
                                for n3 = 0:N
                                    for mu3 = -n3:n3
                                        % Update indices in A and U2 mat
                                        q3 = n3*(n3+1)+mu3;
                                        
                                        tmpSum = 0;
                                        
                                        for l = 0:N
                                            for m = -l:l
                                                % Update indices in A and U2 mat
                                                p = l*(l+1)+m;
                                                
                                                tmpSum = tmpSum + ...
                                                    G(p+1, q1+1, q2+1)*G(q+1, p+1, q3+1);
                                            end
                                        end
                                        A(q+1,q1+1,q2+1,q3+1) = tmpSum;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
    case 5
        A = zeros((N+1)^2, (N+1)^2, (N+1)^2, (N+1)^2, (N+1)^2);
        for n = 0:N
            for mu = -n:n
                % Update indices in A and U mat
                q = n*(n+1)+mu;
                
                for n1 = 0:N
                    if n1>0; fprintf('n=%d n1=%d n=2%d n3=%d n4=%d\n',n,n1,n2,n3,n4); end
                    for mu1 = -n1:n1
                        % Update indices in A and U1 mat
                        q1 = n1*(n1+1)+mu1;
                        
                        for n2 = 0:N
                            for mu2 = -n2:n2
                                % Update indices in A and U2 mat
                                q2 = n2*(n2+1)+mu2;
                                
                                for n3 = 0:N
                                    for mu3 = -n3:n3
                                        % Update indices in A and U2 mat
                                        q3 = n3*(n3+1)+mu3;
                                        
                                        for n4 = 0:N
                                            for mu4 = -n4:n4
                                                % Update indices in A and U2 mat
                                                q4 = n4*(n4+1)+mu4;
                                                
                                                tmpSum = 0;
                                                
                                                for l = 0:N
                                                    for m = -l:l
                                                        % Update indices in A and U2 mat
                                                        p = l*(l+1)+m;
                                                        
                                                        for l1 = 0:N
                                                            for m1 = -l1:l1
                                                                % Update indices in A and U2 mat
                                                                p1 = l1*(l1+1)+m1;
                                                                
                                                                tmpSum = tmpSum + ...
                                                                    G(p1+1, q2+1, q3+1)*G(p+1, p1+1, q1+1)*G(q+1, p+1, q4+1);
                                                            end
                                                        end
                                                    end
                                                end
                                                A(q+1,q1+1,q2+1,q3+1,q4+1) = tmpSum;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
end


function A = real_gaunt_mtx(N1, N2, N)
%REAL_GAUNT_MTX Construct a matrix of Gaunt coefficients
%
% REAL_GAUNT_MTX constructs the (N1+1)^2x(N2+1)^2x(N+1)^2 matrix of Real 
% Gaunt (R-Gaunt) coefficients which represent the instegral of three 
% spherical harmonics such as
% G_R^q_{q',q''} = \int_\Omega Y_{q'}Y_{q''}Y_{q} \mathrm{d}\Omega.
%
% The R-Gaunt coefficients are calculated from the complex ones using the
% equation(26) in https://doi.org/10.1016/S0166-1280(96)90531-X
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Michele Guerreri (m.guerreri@ucl.ac.uk)
%   adapted from:
%   https://github.com/polarch/Spherical-Harmonic-Transform
%   Archontis Politis, 10/10/2013
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute GAUNT coeff from complex set:

A = zeros((N1+1)^2, (N2+1)^2, (N+1)^2);
for n = 0:N
    % Compuute unitary matrix to switch from complex 2 real SH set
    U = Im2Real_SH(n);
    % compute complex conjugate of U
    Ucc = U';
        % initialize counter for U matrix rows
        cq = 0;
        for mu = -n:n
            % Update indices in A and U mat
            q = n*(n+1)+mu;
            cq = cq + 1;
                        
            for n1 = 0:N1
                % Compuute unitary matrix to switch from complex 2 real SH set
                U1 = Im2Real_SH(n1);
                % initialize counter for U1 matrix rows
                cq1 = 0;
                for mu1 = -n1:n1
                    % Update indices in A and U1 mat
                    q1 = n1*(n1+1)+mu1;
                    cq1 = cq1 + 1;
                    
                    for n2 = 0:N2
                        % Compuute unitary matrix to switch from complex 2 real SH set
                        U2 = Im2Real_SH(n2);
                        % initialize counter for U2 matrix rows
                        cq2 = 0;
                        for mu2 = -n2:n2
                            % Update indices in A and U2 mat
                            q2 = n2*(n2+1)+mu2;
                            cq2 = cq2 + 1;
                            
                            if n<abs(n1-n2) || n>n1+n2
                                %A(q1+1, q2+1, q+1) = 0;
                                break
                            else
                                % temporary summation term
                                tmpSum = 0;
                                % initialize counter for U matrix colums
                                cp = 0;
                                for m = -n:n
                                    % Update indices in G and U mat
                                    cp = cp + 1;
                                    % initialize counter for U1 matrix colums
                                    cp1 = 0;
                                    for m1 = -n1:n1
                                        % Update indices in G and U1 mat
                                        cp1 = cp1 + 1;
                                        % initialize counter for U2 matrix colums
                                        cp2 = 0;
                                        for m2 = -n2:n2
                                            % Update indices in G and U2 mat
                                            cp2 = cp2 + 1;
                                            
                                            % compute the complex Gaunt coeff
                                            wigner3jm = w3j(n1, n2, n, m1, m2, -m);
                                            wigner3j0 = w3j(n1, n2, n, 0, 0, 0);
                                            G = (-1)^m * sqrt((2*n1+1)*(2*n2+1)*(2*n+1)/(4*pi)) * wigner3jm * wigner3j0;
                                            
                                            % Update the sum
                                            tmpSum = tmpSum + ...
                                                G*Ucc(cp,cq)*U1(cq1,cp1)*U2(cq2,cp2);                                                
                                        end
                                    end
                                end
                                % Define the real coeff
                                A(q1+1, q2+1, q+1) = real(tmpSum);
                            end
                        end
                    end
                end
            end
        end
end

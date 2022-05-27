function [l1,l2,l3,v1,v2,v3] = DTIfit(E, bvals, bvecs, fitS0)
% 
% Compute the diffusion tensor from diffusion signal, bvalues and bvecs
% 
% 

% define the outputs
n_meas = size(E,1);
l1 = zeros(n_meas, 1);
l2 = zeros(n_meas, 1);
l3 = zeros(n_meas, 1);
v1 = zeros(n_meas, 3);
v2 = zeros(n_meas, 3);
v3 = zeros(n_meas, 3);


% first get the design matrix
X = [ones(1, length(bvals)); ...
    -bvals'.*bvecs(:,1)'.*bvecs(:,1)'; ...
    -2*bvals'.*bvecs(:,1)'.*bvecs(:,2)'; ...
    -2*bvals'.*bvecs(:,1)'.*bvecs(:,3)'; ...
    -bvals'.*bvecs(:,2)'.*bvecs(:,2)'; ...
    -2*bvals'.*bvecs(:,2)'.*bvecs(:,3)'; ...
    -bvals'.*bvecs(:,3)'.*bvecs(:,3)']';

% Get design matrix pseudo-inverse
if fitS0 == 0
  X = X(:,2:end);
end
Xi = pinv(X);

% do the fit
if (fitS0)
	D = Xi*log(E)';
else
	S0 = squeeze(mean(E(:,bvals==0),2));
	E = E./S0;
	D = zeros(7,n_meas);
	D(1,:) = log(S0);
	D(2:end,:) = Xi*log(E)';
end

% Reorganize the output
D = reshape(real(D), [size(D,1) 1 size(D,2)]);

DT = [  D(2,:,:) D(3,:,:) D(4,:,:) ;  ...
        D(3,:,:) D(5,:,:) D(6,:,:) ; ...
        D(4,:,:) D(6,:,:)  D(7,:,:)];

% get eigenvalues and eigenvectors    
for i = 1:n_meas
    [V, L] = eig(DT(:,:,i));
    % get orderd eigenvalues
    [LOrd, LOrd_idx] = sort( diag(L), 'descend');
    % Job done !
    v1(i,:) = V(:,LOrd_idx(1));
    v2(i,:) = V(:,LOrd_idx(2));
    v3(i,:) = V(:,LOrd_idx(3));
    l1(i) = LOrd(1);
    l2(i) = LOrd(2);
    l3(i) = LOrd(3);
end



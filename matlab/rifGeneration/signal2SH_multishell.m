function shCoeff = signal2SH_multishell(signal, protocol, lmax)
% 
% Computes the SH coefficients from the input signal
%
% Usage: 
%   shCoeff = signal2SH_multishell(signal, protocol)
%   shCoeff = signal2SH_multishell(signal, protocol, lmax)
%
% input:
%   signal      Nvox x Nsignal array 
%   protocol    protocol under the NODDI matlab toolbox framework.
%   lmax        should be an array of lmax vaues, one for each 
%               combination of b-vaslues and b-shapes.
% 
% Output:
%   shCoeff       ...
% 
% Author:
%   Michele Guerreri (m.guerreri@ucl.ac.uk)
% 

% How many unique b-values are there?
b = GetB_Values(protocol);
if strcmp(protocol.pulseseq, 'CTE')
    b_nn0 = b;
    b_nn0(protocol.b0_Indices) = [];
    a_nn0 = protocol.alpha;
    a_nn0(protocol.b0_Indices) = [];
    ub = unique( cat(1, b_nn0, a_nn0)',"rows" );
else
    ub = Get_uB_Values(protocol)';
end
n_ub = size(ub,1);

% make sure ub and lmax are consisent
if n_ub ~= size(lmax,1)
    error('Mismatch between the number of lmax input and the b-shells found in this protocol\n');
end

% define the output
shCoeff = cell(n_ub, 1);

% loop over the shells
for ii = 1:n_ub
    % Determine the max order based on the number of directions in that shell
    if strcmp(protocol.pulseseq, 'CTE')
        tmp_b_idx = find( b == ub(ii,1) & protocol.alpha == ub(ii,2) );
    else
        tmp_b_idx = find( b == ub(ii));
    end
    nDir = length(tmp_b_idx);
    % select signal and direction
    S = signal(:,tmp_b_idx);
    dirs = protocol.grad_dirs(tmp_b_idx,:);
    % Select lmax
    l = lmax(ii);
        
    % communicate with the user
    fprintf('--------------------------------------------------------------\n')
    fprintf('Calculating the rotated SH coefficients of b-shell %.0f s/mm^2.\n', ...
        ub(ii)*1e-6)
    if strcmp(protocol.pulseseq, 'CTE')
        fprintf('CTE pulse sequence is used with a b-tensor shape of %.2f.\n', ...
        ub(ii,2))
    end
    fprintf('There are %i directions and %i degrees of freedom.\n', ...
        nDir, (l+1)*(l/2+1))

    % Compute the SH for this shell
    sh_scheme = gen_scheme(dirs, l);
    shCoeff{ii} = amp2SH(S', sh_scheme)';
end

    

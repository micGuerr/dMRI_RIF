function [V1,V2,V3] = getDTbasis(signal, protocol, varargin)
% 
% This function compute the diffusion tensor orthonormal basis from a set
% of diffusion weighted signals.
% 
% Usage: 
%   [V1,V2,V3] = getDTbasis(signal, protocol)
%   [V1,V2,V3] = getDTbasis(signal, protocol, bvals)
%
% input:
%   signal      Nvox x Nsignal matrics 
%   protocol    protocol under the NODDI matlab toolbox framework.
%   bvals       [optional] specify the b-value to use in s/mm2. Default is
%               the closest b-vaslue to 1200 s/mm2
% 
% Output:
%   [V1,V2,V3]  Nvox x 3 indicating the three basis element of the
%               diffusion tensor
% 
% Author:
%   Michele Guerreri (m.guerreri@ucl.ac.uk)
% 

% get the b-values
b = GetB_Values(protocol)'*1e-9;


% if bvals not spec
if isempty(varargin)
    % a b-value which is the closest to 1200 s/mm2. This
    % should ensure that independently from the protocol we use a good b-vale
    % to calulate DTI
    bu = unique(b);
    [~, dti_bIdx ] = min( abs(bu - 1.3) );
    dti_b =  bu(dti_bIdx);
    if strcmp(protocol.pulseseq, 'CTE')
        b_idx_nn0 = find(b == dti_b & protocol.alpha' == 1);
    else
        b_idx_nn0 = find(b == dti_b);
    end
else
    dti_b = varargin{1}*1e-3;
    if strcmp(protocol.pulseseq, 'CTE')
        b_idx_nn0 = find( round(b,3) == round(dti_b,3) &  ...
                          protocol.alpha' == 1);
    else
        b_idx_nn0 = find( round(b,3) == round(dti_b,3) );
    end
end
% add the b=0 measure for dti computation
b_idx = cat(1, protocol.b0_Indices ,b_idx_nn0);

% communiocate with the user
 fprintf('\nFitting DT using b-value %.0f s/mm^2 and number of directions %i...\n', ...
     dti_b*1e3, length(b_idx_nn0) );

% define the entries for DT computations
E = signal(:, b_idx);
BVL = b(b_idx);
BVC = protocol.grad_dirs(b_idx,:);

% now compute the dt
[~,~,~,V1,V2,V3] = DTIfit(E, BVL, BVC, 0);







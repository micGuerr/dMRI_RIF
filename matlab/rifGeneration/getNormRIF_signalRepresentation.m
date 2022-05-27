function rifnorm = getNormRIF_signalRepresentation(signal, protocol, lmax, d)
%
% Compute rotation invariant features from a SH representation of the
% signal given a specific set of lmax for each unique acquisition set
% (b-value and b-shape).
%
% Usage:
% rifnorm = getNormRIF_signalRepresentation(signal, protocol, lmax, d)
% 
% Input:
%   signal      NxM array of signal. N is the number of voxels and M is the
%               number of acquisitions.
%   protocol    struct specifying the acquisition protocol.
%   lmax        Kx1 array specifying the lmax to use in the expansion for
%               each uniqe acquisition type. The unicity is given by the 
%               b-value and b-shape used in the acquisition.
%   d           Kx1 array specifying the power of the invariants.
% 
% Output:
% rifnorm    MxW array of SH coeff. representing the signal
%
% Author
%   Michele Guerreri
%
%

% load the indices for alebrical independence
load('C:\Users\ucacmgu\OneDrive - University College London\NODDIplus-DNN\Utilities\algIndp.mat', 'algIndpRif')

% compute the SH coefficients
SHcoeff = signal2SH_multishell(signal, protocol, lmax);

% get the main fibre direction
V1 = getDTbasis(signal, protocol);

% temporary storage of the rif and their normalization
tmp_rif = cell( numel(SHcoeff),max(d) );
tmp_norm = cell( numel(SHcoeff),max(d) );


for ii = 1:numel(SHcoeff)
    delta = gen_scheme(V1, lmax(ii));
    for jj = 1:d(ii)
        tmp_rif{ii,jj} =  getRIF_fullVec(SHcoeff{ii}, lmax(ii), jj);
        tmp_norm{ii,jj} =  getRIF_fullVec(delta.sh, lmax(ii), jj);
    end
end


% declear the output
rif = [];
normalization = [];
idx = [];
% loop over the shells to concatenate SHs coeff.
for ii = 1:numel(SHcoeff)
    for jj = 1:d(ii)
        rif = cat( 2, rif, tmp_rif{ii,jj});
        normalization = cat( 2, normalization, tmp_norm{ii,jj});
    end
    idx = cat( 1, idx, algIndpRif{lmax(ii)/2,d(ii)});
end
rif = rif(:,idx==1);
normalization = normalization(:,idx==1);
rifnorm = rif./normalization;


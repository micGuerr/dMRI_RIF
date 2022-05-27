function rif = getRIF_signalRepresentation(signal, protocol, lmax, d)
%
% Compute rotation invariant features from a SH representation of the
% signal given a specific set of lmax for each unique acquisition set
% (b-value and b-shape).
%
% Usage:
% rif_signal = getRIF_signalRepresentation(signal, protocol, lmax, d)
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
% rif_signal    MxW array of SH coeff. representing the signal
%
% Author
%   Michele Guerreri
%
%


% compute the SH coefficients
SHcoeff = signal2SH_multishell(signal, protocol, lmax);

% temporary storage of the rif
tmp_rif = cell( numel(SHcoeff),1 );

for ii = 1:numel(SHcoeff)
    tmp_rif{ii} =  getRIF(SHcoeff{ii}, lmax(ii), d(ii));
end


% declear the output
rif = [];
% loop over the shells to concatenate SHs coeff.
for ii = 1:numel(SHcoeff)
    rif = cat( 2, rif, tmp_rif{ii});
end



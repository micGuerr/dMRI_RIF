function protocol = genNODDI_protocol(ubval, n_bdir, varargin)
% 
% Generate a protocol structure compatible with NODDI Matlab toolbox from
% the specified set of b-values and number of directions. The shape o the
% btensor can lso be input
%
% Usage: 
%   protocol = genNODDI_protocol(ubval, n_bdir)
%   protocol = genNODDI_protocol(ubval, n_bdir, ubshape)
%
% input:
%   ubval       1xN array of unique b-values in s/mm2.  
%               For example [0 1000 2000]
%   n_bdir      Number of direction to be used for each b-values.
%               Must be consistent with bval.
%               For example [5 30 60]
%   ubshape     OPTIONAL. shape of the b-tensor
%               For example [1 1 0.5]
% 
% Output:
%   protocol    protocol structure compatible with NODDI Matlab toolbox
% 
% Requirements:
%   CSD matlab toolbox (https://github.com/jdtournier/csd)
%
% Author:
%   Michele Guerreri (m.guerreri@ucl.ac.uk)
% 

if nargin > 2
    ubshape = varargin{1};
else
    ubshape = [];
end

% if bshape file is provided, pulse seq is Tensor-values. Else it is PGSE
if ~isempty(ubshape)
    % if the bshape file is provided, use "Cylindrical Tensor Encoding" pulse seq
    protocol.pulseseq = 'CTE';
else
    protocol.pulseseq = 'PGSE';
end
protocol.schemetype = 'multishellfixedG';
protocol.teststrategy = 'fixed';

%% make sure all inputs have same length

if size(ubval,2) ~= size(n_bdir,2)
    error('mismatch between number of b-values and directions per b-shell.')
end
if ~isempty(ubshape)
    if size(ubval,2) ~= size(ubshape,2)
        error('mismatch between number of b-values and type of shapes.')
    end
end

%% get bvalues
n_shells = size(ubval,2);

bval = [];
for i = 1:n_shells
    bval = cat(1, bval, repmat(ubval(i), [n_bdir(i) 1]) );
end

% if the b0 threshold should be different from 0 change this
b0threshold = 0;

% set total number of measurements
protocol.totalmeas = length(bval);

% set the b=0 indices
protocol.b0_Indices = find(bval<=b0threshold);
protocol.numZeros = length(protocol.b0_Indices);

% find the unique non-zero b-values
B = unique(bval(bval>b0threshold));
% set the number of shells
protocol.M = length(B);
for i=1:length(B)
    protocol.N(i) = length(find(bval==B(i)));
end

% maximum b-value in the s/mm^2 unit
maxB = max(B);

% set maximum G = 40 mT/m
Gmax = 0.04;

% set smalldel and delta and G
GAMMA = 2.675987E8;
tmp = nthroot(3*maxB*10^6/(2*GAMMA^2*Gmax^2),3);
for i=1:length(B)
    protocol.udelta(i) = tmp;
    protocol.usmalldel(i) = tmp;
    protocol.uG(i) = sqrt(B(i)/maxB)*Gmax;        
end

protocol.delta = zeros(size(bval))';
protocol.smalldel = zeros(size(bval))';
protocol.G = zeros(size(bval))';

for i=1:length(B)
    tmp = find(bval==B(i));
    for j=1:length(tmp)
        protocol.delta(tmp(j)) = protocol.udelta(i);
        protocol.smalldel(tmp(j)) = protocol.usmalldel(i);
        protocol.G(tmp(j)) = protocol.uG(i);
    end
end

%% get bvecs

bvec = [];
for i = 1:n_shells
    if ubval(i) == 0
        bvec = cat(1, bvec, zeros(n_bdir(i), 3) );
    else
        bvec = cat(1, bvec, equidistribute(n_bdir(i)) );
    end
end
protocol.grad_dirs = bvec;

% some systems set vector to zeros for b=0
% the codes below try to account for this
if isempty(protocol.b0_Indices)
    for i=1:protocol.totalmeas
        if norm(protocol.grad_dirs(i,:)) == 0
            protocol.G(i) = 0.0;
            protocol.b0_Indices = [protocol.b0_Indices i];
        end
    end
    protocol.numZeros = length(protocol.b0_Indices);
end
    
% make the gradient directions for b=0's [1 0 0]
for i=1:length(protocol.b0_Indices)
    protocol.grad_dirs(protocol.b0_Indices(i),:) = [1 0 0];
end

% make sure the gradient directions are unit vectors
for i=1:protocol.totalmeas
    protocol.grad_dirs(i,:) = protocol.grad_dirs(i,:)/norm(protocol.grad_dirs(i,:));
end


%% load bshape if needed
if ~isempty(ubshape)
    
    
    
    bshape = [];
    for i = 1:n_shells
        bshape = cat(1, bshape, repmat(ubshape(i), [n_bdir(i) 1]) );
    end



    % find the unique b-shapes
    A = unique(bshape);
    
    % set the number of shapes
    protocol.Ms = length(A);
    for i=1:length(A)
        protocol.Ns(i) = length(find(bshape==A(i)));
        protocol.ushape(i) = A(i);
    end
    % Store the bshapes per acquisition
    protocol.alpha = bshape';
end







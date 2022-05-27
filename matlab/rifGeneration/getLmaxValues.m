function lmax = getLmaxValues(protocol, rep_type)
%
% Finds the optimal lmax to use for SH expansion of the diffusion signal
% given an acquisition protocol and a specific degree of flexibility.
%
% Usage:
% lmax = getLmaxValues(protocol, rep_type)
% 
% Input:
%   protocol    is the acquisition protocol in NODDI-toolbox format.
%   rep_type    is a string specifying the flexibility required in the SH
%               estimation.
%               we consider four types of SH representations:
%               1. highly overdetermined (number of SH coeff << number of direction).
%               2. moderately overdetermined (number of SH coeff < number of direction).
%               3. littlest overdetermined (number of SH coeff <~ number of direction).
%               4. underdetermined (number of SH coeff > number of direction).
% 
% Output:
% lmax          array of lmax values, one for each unique acquisition type. The
%               unicity is given by the b-value and b-shape used in the
%               acquisition.
%
% Author
%   Michele Guerreri
%

% A list of rhe degrees of freedom corresponding to the first 2N degrees
SHdof = ((1:30).*2+1).*((1:30)+1);

% how many b-values and b-shapes are there?
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
    ub = cat(2,  ub, ones(length(ub),1));
end
n_ub = size(ub,1);

% let's first define output
lmax = nan(n_ub,1);


% loop over the b-shells
for bb = 1:n_ub
    % Determine the number of directions in that shell
    if strcmp(protocol.pulseseq, 'CTE')
        tmp_b_idx = find( b == ub(bb,1) & protocol.alpha == ub(bb,2) );
    else
        tmp_b_idx = find( b == ub(bb));
    end
    nDir = length(tmp_b_idx);
    % check that there are at least 6 directions
    if nDir < 6
        error(' in order to have an overdetermined SH representation at least 6 directions must be input\n');
    end
    
    % Find the index for which the number of directions of the shell is just
    % under the SH degrees of freedom
    nDir_SHdof_ratioLog = log( nDir./SHdof );
    [ ~, N ] = min(nDir_SHdof_ratioLog(nDir_SHdof_ratioLog >= 0 ));
    
    
    % use SH degree > 2 only if the signal is LTE
    if ub(bb,2) == 1
        % swhitch between the representation types and assigne the
        % opporpriate l max
        switch rep_type
            case 'hover'
                lmax(bb) = 2;
            case 'mover'
                if N > 1
                    lmax(bb) = (N-1)*2;
                else
                    lmax(bb) = N*2;
                end
            case 'lover'
                lmax(bb) = N*2;
            case 'under'
                lmax(bb) = (N+1)*2;
        end
    else
        lmax(bb) = 2;
    end
end








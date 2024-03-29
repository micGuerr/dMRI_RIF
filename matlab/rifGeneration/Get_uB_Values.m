function b = Get_uB_Values(protocol)
% Computes the b value for each measurement in the protocol
%
% b = Get_uB_Values(protocol)
% returns an array of unique b values, one for each unique measurement defined
% in the list of measurements in the protocol.
%
% author: 
%   Daniel C Alexander (d.alexander@ucl.ac.uk)
%   Michele Guerreri (m.guerreri@ucl.ac.uk)
%

GAMMA = 2.675987E8;
if(strcmp(protocol.pulseseq, 'PGSE') || strcmp(protocol.pulseseq, 'STEAM') ...
        || strcmp(protocol.pulseseq, 'CTE'))
    modQ = GAMMA*protocol.usmalldel.*protocol.uG;
    diffTime = protocol.udelta - protocol.usmalldel/3;
    b = diffTime.*modQ.^2;

elseif(strcmp(protocol.pulseseq, 'DSE'))
    b = getB_ValuesDSE(protocol.G, protocol.delta1, protocol.delta2, protocol.delta3, protocol.t1, protocol.t2, protocol.t3);

elseif(strcmp(protocol.pulseseq, 'OGSE'))
    b = GetB_ValuesOGSE(protocol.G, protocol.delta, protocol.smalldel, protocol.omega);
end




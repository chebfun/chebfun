function f = divergence(F)
%DIVERGENCE   Divergence of a CHEBFUN3V object.
%   DIVERGENCE(F) returns divergence of the CHEBFUN3V object F as a 
%   CHEBFUN3. If F = U i + V j + W k, then divergence(F) = U_x + V_y + W_z.
%
% See also CHEBFUN3V/DIV.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) )
    f = chebfun3();
    return
end


if ( F.nComponents == 3 )
    Fc = F.components; 
    diff1 = diff(Fc{1}, 1, 1);
    diff2 = diff(Fc{2}, 1, 2);
    diff3 = diff(Fc{3}, 1, 3);
    f = diff1 + diff2 + diff3;
else                      
     error('CHEBFUN:CHEBFUN3V:divergence:notSupported', ...
        'three inputs are needed.')
end

end
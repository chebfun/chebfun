function F = curl(F)
%CURL   Curl of a CHEBFUN3V object.
%   G = CURL(F) returns the CHEBFUN3V of the curl of F. If F is a CHEBFUN3V 
%   with three components then it returns the CHEBFUN3V representing the 3D
%   curl operation.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) )
    F = chebfun3v;
    return
end

Fc = F.components; 

if ( F.nComponents == 3 )   % 3D curl
   F = [diff(Fc{3}, 1, 2) - diff(Fc{2}, 1, 3);
        diff(Fc{1}, 1, 3) - diff(Fc{3}, 1, 1);
        diff(Fc{2}, 1, 1) - diff(Fc{1}, 1, 2)];
else                      
     error('CHEBFUN:CHEBFUN3V:curl:notSupported', ...
        'Defined only for CHEBFUN3V objects with three components.')
end

end
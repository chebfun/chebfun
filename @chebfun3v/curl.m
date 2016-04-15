function F = curl(F)
%CURL  curl of a CHEBFUN3V
%   S = CURL(F) returns the CHEBFUN3 of the curl of F. If F is a CHEBFUN3V with
%   with three components then it returns the CHEBFUN3V representing the 3D curl operation.

% Empty check: 
if ( isempty( F ) )
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
        'three inputs are needed.')
end

end
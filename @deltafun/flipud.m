function f = flipud(f)
%FLIPUD   Flip/reverse a DELTAFUN object.
%   G = FLIPUD(F) returns G such that G(x) = F(-x) for all x in [-1,1].

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org for Chebfun information.

if ( isempty(f) )
    return
end

% Flip the funPart:
if ( ~isempty(f.funPart) )
    f.funPart = flipud(f.funPart);
    if ( ~isempty(f.deltaLoc) )
        % Map deltaMag to [-1, 1]:
        inverseMap = f.funPart.mapping.inv;
        loc = inverseMap(f.deltaLoc);
        % Location is a vector, so flipud translates into fliplr:
        loc = fliplr(-loc);
        f.deltaMag = fliplr(f.deltaMag);
        % Map back the locations:
        forwardMap = f.funPart.mapping.for;
        f.deltaLoc = forwardMap(loc);
    end
else
    if ( ~isempty(f.deltaLoc) )
        % No map given, reflect about the origin:
        f.deltaLoc = fliplr(-f.deltaLoc);
        f.deltaMag = fliplr(f.deltaMag);
    end
end       
end

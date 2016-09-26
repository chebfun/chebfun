function f = flipud(f)
%FLIPUD   Flip/reverse a DELTAFUN object.
%   G = FLIPUD(F) returns G such that G(x) = F(-x) for all x in [-1,1].

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    return
end

% Flip the funPart:
f.funPart = flipud(f.funPart);
if ( ~isempty(f.deltaLoc) )
    % Map deltaMag to [-1, 1]:
    inverseMap = f.funPart.mapping.Inv;
    loc = inverseMap(f.deltaLoc);
    % Location is a vector, so flipud translates into fliplr:
    loc = fliplr(-loc);
    f.deltaMag = fliplr(f.deltaMag);
    % Map back the locations:
    forwardMap = f.funPart.mapping.For;
    f.deltaLoc = forwardMap(loc);
end

end

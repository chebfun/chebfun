function [g, jumpVals, locations] = cumsum(f)
%CUMSUM   Indefinite integral of a DELTAFUN.
%   CUMSUM(F) is the indefinite integral of the DELTAFUN F. LOCATIONS is a
%   vector which indicates the locations of the delta functions only (not their
%   derivatives). JUMPVALS is the vector of (signed) magnitude of jumps, that 
%   should be introduced at these locations.
% See also SUM

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize output:

jumpVals = [];
locations = [];
g = deltafun;

% Trivial case:
if ( isempty(f) )
    g.funPart = fun.constructor(0);
    return;
end

deltaMag = f.impulses;
deltaLoc = f.location;

if ( isempty(deltaLoc) || isempty(deltaMag) )
    g = cumsum(f.funPart);
    jumpVals = [];
    locations = [];
else
    if ( size(deltaMag, 1) > 1 )
        g.impulses = deltaMag(2:end, :);
        g.location = deltaLoc;
    end
    g = simplify(g);
    
    % Integrate the funPart:
    if ( isempty(f.funPart) )
        g.funPart = fun.constructor(0);
    else                
        g.funPart = cumsum(f.funPart);
    end
    
    pref = chebpref();
    deltaTol = pref.deltaPrefs.deltaTol;
    
    jumpVals = deltaMag(1, :);
    idx = abs(jumpVals) > deltaTol; 
    jumpVals = jumpVals(idx);
    locations = deltaLoc(idx);            
end

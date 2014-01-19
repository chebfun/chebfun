function [g, jumpVals, locations] = cumsum(f)
%CUMSUM   Indefinite integral of a DELTAFUN.
%   CUMSUM(F) is the indefinite integral of the DELTAFUN F.
%

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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
    
    % Get the domain:
    if ( isempty(f.funPart) )
        g.funPart = fun.constructor(0);
    else                
        g.funPart = cumsum(f.funPart);
    end
    
    tol = deltafun.pref.deltafun.deltaTol;
    
    jumpVals = deltaMag(1, :);
    idx = abs(jumpVals) > tol; 
    jumpVals = jumpVals(idx);
    locations = deltaLoc(idx);            
end


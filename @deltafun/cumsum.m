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
    % Return an empty deltafun in this case:
    return;
end

deltaMag = f.deltaMag;
deltaLoc = f.location;

if ( isempty(deltaLoc) || isempty(deltaMag) )
    g = cumsum(f.funPart);    
else    
    g = {};
    pref = chebpref();
    deltaTol = pref.deltaPrefs.deltaTol;
    idx = abs(deltaMag(1,:)) >= deltaTol;
    newBreaks = deltaLoc(idx);
    jumpVals = deltaMag(1, idx);
    breakPts = union(f.funPart.domain, newBreaks);
    funParts = restrict(f.funPart, breakPts);
    for k = 1:numel(funParts)
        fk = cumsum(funParts{k}) + (k-1);
        dk = fk.domain;
        idx = deltaLoc >= dk(1) & deltaLoc <= dk(2);
        lk = deltaLoc(idx);
        mk = deltaMag(2:end, idx);
        if ( isempty(idx) || all(mk == 0) )
            g{k} = fk;
        else
            g{k} = deltafun(fk, dk, lk);
        end
        
    end
    
    % Integrate the funPart:
    if ( isempty(f.funPart) )
        g.funPart = fun.constructor(0);
    else                
        g.funPart = cumsum(f.funPart);
    end
    
    % Get the tolerance:
    pref = chebpref();
    deltaTol = pref.deltaPrefs.deltaTol;
    
    jumpVals = deltaMag(1, :);
    idx = abs(jumpVals) > deltaTol; 
    jumpVals = jumpVals(idx);
    locations = deltaLoc(idx);            
end

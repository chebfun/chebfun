function g = cumsum(f)
%CUMSUM   Indefinite integral of a DELTAFUN.
%   CUMSUM(F) is the indefinite integral of the DELTAFUN F.
%

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If f has no singularity at any endpoint, then just integrate its smooth part
% and return.

deltaMag = f.impulses;
deltaLoc = f.locations;
if ( isempty(deltaMag) || isemtpy(deltaLoc) )
    g = cumsum(f.funPart);
else
    g = deltafun;
    if ( size(deltaMag, 1) > 1 )
        g.impulses = deltaMag(2:end, :);
        g.locations = deltaLoc;
    end
    
    % Get the domain:
    dom = f.funPart.domain;
    tol = deltafun.pref.deltafun.deltaTol;
 
    F = cumsum(f.funPart);
    
    for j = 1:length(deltaLoc)
        if ( abs(deltaMag(1, j)) > tol )
        end
    end
            
end
g = simplify(g);


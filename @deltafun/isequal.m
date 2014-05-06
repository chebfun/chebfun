function out = isequal(f, g)
%ISEQUAL   Test if DELTAFUN objects F and G are equal.
%   ISEQUAL(F, G) returns TRUE if the DELTAFUN objects F and G have the same
%   underlying FUNPART and the same delta functions.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Get the tolerances:
pref = chebfunpref();
proximityTol = pref.deltaPrefs.proximityTol;
deltaTol = pref.deltaPrefs.deltaTol;

% Assume true be default:
out = 1;

% Trivial cases:
if ( isempty(f) && isempty(g) )
    out = 1;
    return
elseif ( xor(isempty(f), isempty(g)) )
    out = 0;
    return
end

% Non-trivial cases
if ( ~isequal(f.funPart, g.funPart) )
    out = 0;
    return
end

% Simplify f and g:
f = simplifyDeltas(f);
g = simplifyDeltas(g);

% Test sizes and values of deltaLoc and deltaMag fields:
if ( any( size(f.deltaLoc) ~= size(g.deltaLoc) ) )
    out = 0;
    return
end
if ( any(abs((f.deltaLoc - g.deltaLoc)) > proximityTol) )
    out = 0;
    return
end
if ( any(size(f.deltaMag) ~= size(g.deltaMag)) )
    out = 0;
    return;
end
if ( any(abs(f.deltaMag - g.deltaMag) > deltaTol) )
    out = 0;
    return
end

end

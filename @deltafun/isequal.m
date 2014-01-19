function out = isequal(f, g)
%ISEQUAL   Test if DELTAFUN objects F and G are equal.
%   ISEQUAL(F, G) returns TRUE if the DELTAFUN objects F and G have the same
%   underlyig FUNPART and the same delta functions.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.


pTol = deltafun.pref.deltafun.proximityTol;
dTol = deltafun.pref.deltafun.deltaTol;

% Assume ture be default:
out = 1;

% Trivial cases:
if ( isempty(f) && isempty(g) )
    out = 1;
    return
end

if ( xor(isempty(f), isempty(g)) )
    out = 0;
    return
end

% Non-trivial cases
if ( ~iszero(f.funPart - g.funPart) )
    out = 0;
    return
end


f = simplify(f);
g = simplify(g);

if ( any( size(f.location) ~= size(g.location) ) )
    out = 0;
    return
end

if ( any(abs((f.location - g.location)) > pTol) )
    out = 0;
    return
end

if ( any(size(f.impulses) ~= size(g.impulses)) )
    out = 0;
    return;
end

if ( any(abs(f.impulses - g.impulses) > dTol) )
    out = 0;
    return
end


end

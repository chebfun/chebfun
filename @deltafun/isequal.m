function out = isequal(f, g)
%ISEQUAL   Test if DELTAFUN objects F and G are equal.
%   ISEQUAL(F, G) returns TRUE if the DELTAFUN objects F and G have the same
%   underlyig FUNPART and the same delta functions.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.


pTol = deltafun.pref.deltafun.proximityTol;
dTol = deltafun.pref.deltafun.deltaTol;


if ( f.funPart ~= g.funPart )
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

if ( any( size(f.magnitude) ~= size(g.magnitude) ) )
    out = 0;
    return;
end

if ( any( size(f.magnitude) ~= size(g.magnitude) ) )
    out = 0;
    return
end


end

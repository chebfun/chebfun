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

if ( any( size(f.delta.location) ~= size(g.delta.location) ) )
    out = 0;
    return
end

if ( any(abs((f.delta.location - g.delta.location)) > pTol) )
    out = 0;
    return
end

if ( any( size(f.delta.magnitude) ~= size(g.delta.magnitude) ) )
    out = 0;
    return;
end


fdeltaMag = f.delta.magnitude;
gdeltaMag = g.delta.magnitude;

if ( any( size(fdeltaMag) ~= size(gdeltaMag) ) )
    out = 0;
    return
end


end

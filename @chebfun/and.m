function h = and(f, g)
%&   CHEBFUN Logical AND
%   F & G performs a logical AND of two CHEBFUN objects F and G and returns a
%   CHEBFUN containing elements set to either logical 1 (TRUE) or logical 0
%   (FALSE). An element of the output CHEBFUN is set to 0 if both input CHEBFUN
%   objects contains a non-zero element at that point, otherwise it is set to 0.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Check for emptiness:
if ( isempty(f) )
    h = f;
    return
elseif ( isempty(g) )
    h = g;
    return
end

% Check the orientation:
isTransposedF = f.isTransposed;
isTransposedG = g.isTransposed;
if ( isTransposedF ~= isTransposedG )
    error('CHEBFUN:or:trans', 'Matrix dimensions must agree.');
elseif ( isTransposedF )
    f = f.';
    g = g.';
end

% Check the domains:
if ( checkDomain )
    error('CHEBFUN:and:doms', 'Inconsistent domains.');
end

% Compute AND() using FLOOR() and ANY():
h = floor(.5*(any(f) + any(g)));

% Transpose result back to row CHEBFUN if required:
if ( isTransposedF )
    h = h.';
end

end
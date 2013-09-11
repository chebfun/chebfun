function h = or(f, g)
%|   CHEBFUN logical OR.
%   F | G performs a logical OR of the CHEBFUN objects F and G and returns a
%   CHEBFUN containing elements set to either logical 1 (TRUE) or logical 0
%   (FALSE).  An element of the output CHEBFUN is set to 1 if either input
%   CHEBFUN contains a non-zero element at that point, otherwise it is set to 0.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for emptiness:
if ( isempty(f) )
    h = f;
    return
elseif ( isempty(g) )
    h = g;
    return
end

% Check the domains:
if ( checkDomain(f, g) )
    error('CHEBFUN:and:doms', 'Inconsistent domains.');
end

% Check the orientation
isTransposedF = f.isTransposed;
isTransposedG = g.isTransposed;
if ( isTransposedF ~= isTransposedG )
    error('CHEBFUN:or:trans', 'Matrix dimensions must agree.');
elseif ( isTransposedF )
    f = f.';
    g = g.';
end

h = logical(logical(f) + logical(g));

% Transpose back to row CHEBFUN:
if ( isTransposedF )
    h = transpose(h);
end

end
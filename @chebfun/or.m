function h = or(f, g)
%|   Chebfun logical OR.
%
% F | G performs a logical OR of chebfun F and G and returns a chebfun
% containing elements set to either logical 1 (TRUE) or logical 0 (FALSE).  An
% element of the output chebfun is set to 1 if either input chebfun contains a
% non-zero element at that point, otherwise it is set to 0.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Check for emptiness
if ( isempty(f) )
    h = f;
    return
elseif ( isempty(g) )
    h = g;
    return
end

% Check the domains
dF = f.domain([1 end]);
dG = g.domain([1 end]);
if ( any(dF ~= dG) )
    error('CHEBFUN:and:doms', 'Inconsistent domains.');
end

% Check the orientation
tF = f.isTransposed;
tG = g.isTransposed;
if ( tF ~= tG )
    error('CHEBFUN:or:trans', 'Matrix dimensions must agree.');
elseif ( tF )
    f = transpose(f);
    g = transpose(g);
end

h = logical(logical(f) + logical(g));

% Transpose back to row chebfun:
if ( tF )
    h = transpose(h);
end

end
function f = and(f, g)
%&   CHEBTECH logical AND.
%   F & G performs a logical AND of two CHEBTECH objects F and G and returns a
%   CHEBTECH containing elements set to either logical 1 (TRUE) or logical 0
%   (FALSE). An element of the output CHEBTECH is set to 1 if both input
%   CHEBTECH objects have a non-zero element at that point, otherwise it is set
%   to 0.  F and G must either be identically zero or have roots in their
%   domains.  If this is not the case, garbage is returned with no warning.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

arbitraryPoint = 0.1273881594;
f.coeffs = feval(f, arbitraryPoint) & feval(g, arbitraryPoint);

end

function g = simplify(f, tol)
%SIMPLIFY  Simplifies the smoothPart of a SINGFUN object.
%  G = SIMPLIFY(F) returns a SINGFUN which has a simplified smoothPart by
%  calling SIMPLIFY in SMOOTHFUN.
%
%  G = SIMPLIFY(F, TOL) does the same thing but uses TOL supplied by the user
%  and pass it to SIMPLIFY in SMOOTHFUN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 2 )
    g.smoothPart = simplify(f.smoothPart, tol);
elseif ( nargin == 1 )
    g.smoothPart = simplify(f.smoothPart);
end

end
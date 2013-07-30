function g = ceil(f)
%CEIL   Pointwise ceiling of a chebfun.
%   G = CEIL(F) returns the chebfun G such that G(x) = CEIL(F(x)) for each x in
%   the domain of F.
%
% See also CHEBFUN/FLOOR, CHEBFUN/ROUND, CHEBFUN/FIX, CEIL.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Deal with unbounded functions:
if ( ~isfinite(f) )
    error('CHEBFUN:ceil:inf', ...
        'Ceil is not defined for functions which diverge to infinity.');
end

g = -floor(-f);

end

function g = acot(f, pref)
%ACOT   Inverse cotangent of a chebfun.
%
% See also COT, ACOTD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @acot, pref);

end
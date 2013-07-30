function g = acosh(f, pref)
%ACOSH   Inverse hypoerbolic cosine of a chebfun.
%
% See also COSH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @acosh, pref);

end
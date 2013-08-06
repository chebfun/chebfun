function g = acotd(f, pref)
%ACOTD   Inverse cotangent of a chebfun, result in degrees.
%
% See also COTD, ACOT.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @acotd, pref);

end
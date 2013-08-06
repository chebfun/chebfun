function g = asinh(f, pref)
%ASINH   Inverse hyperbolic sine of a chebfun.
%
% See also SINH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @asinh, pref);

end
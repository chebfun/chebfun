function g = tanh(f, pref)
%TANH   Hyperbolic tangent of a chebfun.
%
% See also ATAN, TAND.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @tanh, pref);

end
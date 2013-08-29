function g = sinh(f, pref)
%SINH   Hyperbolic sine of a CHEBFUN.
%   SINH(F) computes the hyperbolic sine of the CHEBFUN F.
%
%   SINH(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also ASINH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @sinh, pref);

end

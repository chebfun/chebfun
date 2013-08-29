function g = coth(f, pref)
%COTH   Hyperbolic cotangent of a CHEBFUN.
%   COTH(F) computes the hyperbolic cotangent of the CHEBFUN F.
%
%   COTH(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also ACOTH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @coth, pref);

end
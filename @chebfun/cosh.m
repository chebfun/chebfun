function g = cosh(f, pref)
%COSH   Hyperbolic cosine of a CHEBFUN.
%   COSH(F) computes the hyperbolic cosine of the CHEBFUN F.
%
%   COSH(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also COS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @cosh, pref);

end

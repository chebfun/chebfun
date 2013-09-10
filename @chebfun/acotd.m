function g = acotd(f, pref)
%ACOTD   Inverse cotangent of a CHEBFUN, result in degrees.
%   ACOTD(F) computes the inverse cotangent (in degrees) of the CHEBFUN F.
%
%   ACOTD(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
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

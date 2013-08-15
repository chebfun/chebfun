function g = asind(f, pref)
%ASIND   Inverse sine of a CHEBFUN, result in degrees.
%   ASIND(F) computes the inverse sine (in degrees) of the CHEBFUN F.
%
%   ASIND(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also SIND, ASIN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @asind, pref);

end

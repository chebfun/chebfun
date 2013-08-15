function g = asin(f, pref)
%ASIN   Inverse sine of a CHEBFUN.
%   ASIN(F) computes the inverse sine of the CHEBFUN F.
%
%   ASIN(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also SIN, ASIND.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @asin, pref);

end

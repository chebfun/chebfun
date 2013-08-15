function g = acosh(f, pref)
%ACOSH   Inverse hypoerbolic cosine of a CHEBFUN.
%   ACOSH(F) computes the inverse hypoerbolic cosine of the CHEBFUN F.
%
%   ACOSH(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also COSH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @acosh, pref);

end

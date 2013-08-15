function g = acoth(f, pref)
%ACOTH   Inverse hyperbolic cotangent of a CHEBFUN.
%   ACOTH(F) computes the inverse hyperbolic cotangent of the CHEBFUN F.
%
%   ACOTH(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also COTH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @acoth, pref);

end

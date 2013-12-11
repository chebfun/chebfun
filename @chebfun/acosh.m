function F = acosh(F, pref)
%ACOSH   Inverse hyperbolic cosine of a CHEBFUN.
%   ACOSH(F) computes the inverse hypoerbolic cosine of the CHEBFUN F.
%
%   ACOSH(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also COSH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Loop over the columns of F:
% Call the compose method:
F = compose(F, @acosh, pref);

end

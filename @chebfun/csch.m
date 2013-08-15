function g = csch(f, pref)
%CSCH   Hyperbolic cosecant of a CHEBFUN.
%   CSCH(F) computes the hyperbolic cosecant of the CHEBFUN F.
%
%   CSCH(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also ACSCH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @csch, pref);

end
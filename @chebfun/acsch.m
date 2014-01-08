function g = acsch(f, pref)
%ACSCH   Inverse hyperbolic cosecant of a CHEBFUN.
%   ACSCH(F) computes the inverse hyperbolic cosecant of the CHEBFUN F.
%
%   ACSCH(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also CSCH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Call the compose method:
g = compose(f, @acsch, pref);

end

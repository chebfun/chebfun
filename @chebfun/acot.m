function g = acot(f, pref)
%ACOT   Inverse cotangent of a CHEBFUN.
%   ACOT(F) computes the inverse cotangent of the CHEBFUN F.
%
%   ACOT(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also COT, ACOTD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Call the compose method:
g = compose(f, @acot, pref);

end

function F = coth(F, pref)
%COTH   Hyperbolic cotangent of a CHEBFUN.
%   COTH(F) computes the hyperbolic cotangent of the CHEBFUN F.
%
%   COTH(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also ACOTH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Call the compose method:
F = compose(F, @coth, pref);

end
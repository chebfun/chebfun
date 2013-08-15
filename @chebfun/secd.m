function g = secd(f, pref)
%SECD   Secant of a CHEBFUN, result in degrees.
%   SECD(F) computes the secant (in degrees) of the CHEBFUN F.
%
%   SECD(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also ASECD, SEC.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @secd, pref);

end

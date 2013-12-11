function F = asecd(F, pref)
%ASECD   Inverse secant of a CHEBFUN, result in degrees.
%   ASECD(F) computes the inverse secant (in degrees) of the CHEBFUN F.
%
%   ASECD(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also SECD, ASEC.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Call the compose method:
F = compose(F, @asecd, pref);

end

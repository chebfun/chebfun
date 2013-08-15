function g = asec(f, pref)
%ASEC   Inverse secant of a CHEBFUN.
%   ASEC(F) computes the inverse secant of the CHEBFUN F.
%
%   ASEC(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also SEC, ASECD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @asec, pref);

end

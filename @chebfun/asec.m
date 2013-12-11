function F = asec(F, pref)
%ASEC   Inverse secant of a CHEBFUN.
%   ASEC(F) computes the inverse secant of the CHEBFUN F.
%
%   ASEC(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also SEC, ASECD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Call the compose method:
F = compose(F, @asec, pref);

end

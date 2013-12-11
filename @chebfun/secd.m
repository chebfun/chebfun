function F = secd(F, pref)
%SECD   Secant of a CHEBFUN, result in degrees.
%   SECD(F) computes the secant (in degrees) of the CHEBFUN F.
%
%   SECD(F, PREF) does the same but uses the CHEBPREF object PREF when computing
%   the composition.
%
% See also ASECD, SEC.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Call the compose method:
F = compose(F, @secd, pref);

end

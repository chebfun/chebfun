function g = sec(f, pref)
%SEC   Secant of a CHEBFUN.
%   SEC(F) computes the secant of the CHEBFUN F.
%
%   SEC(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also ASEC, SECD.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Call the compose method:
g = compose(f, @sec, pref);

end
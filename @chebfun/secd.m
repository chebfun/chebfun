function g = secd(f, pref)
%SECD   Secant of a chebfun, result in degrees.
%
% See also ASECD, SEC.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @secd, pref);

end
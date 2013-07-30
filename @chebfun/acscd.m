function g = acscd(f, pref)
%ACSCD   Inverse cosecant of a chebfun, result in degrees.
%
% See also CSCD, ACSC.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @acscd, pref);

end
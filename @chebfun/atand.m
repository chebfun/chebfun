function g = atand(f, pref)
%ATAN   Inverse tangent of a CHEBFUN, result in degrees.
%   ATAN(F) computes the inverse tangent (in degrees) of the CHEBFUN F.
%
%   ATAN(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also TAND, ATAN2D, ATAN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @atand, pref);

end

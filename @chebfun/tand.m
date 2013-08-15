function g = tand(f, pref)
%TAND   Tangent of a CHEBFUN, result in degrees.
%   TAND(F) computes the tangent (in degrees) of the CHEBFUN F.
%
%   TAND(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also ATAND, TAN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @tand, pref);

end

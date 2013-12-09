function g = atan(f, pref)
%ATAN   Inverse tangent of a CHEBFUN.
%   ATAN(F) computes the inverse tangent of the CHEBFUN F.
%
%   ATAN(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also TAN, ATAND.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% [TODO]: atan2.m, atan2d.m

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @atan, pref);

end

function F = asinh(F, pref)
%ASINH   Inverse hyperbolic sine of a CHEBFUN.
%   ASINH(F) computes the inverse hyperbolic sine of the CHEBFUN F.
%
%   ASINH(F, PREF) does the same but uses the CHEBPREF object PREF when
%   computing the composition.
%
% See also SINH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Call the compose method:
F = compose(F, @asinh, pref);

end

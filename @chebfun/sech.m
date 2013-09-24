function g = sech(f, pref)
%SECH   Hyperbolic secant of a CHEBFUN.
%   SECH(F) computes the hyperbolic secant of the CHEBFUN F.
%
%   SECH(F, PREF) does the same but uses the preference structure PREF when
%   computing the composition.
%
% See also ASECH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @sech, pref);

end
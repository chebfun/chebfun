function F = sinh(F, pref)
%SINH   Hyperbolic sine of a CHEBFUN.
%   SINH(F) computes the hyperbolic sine of the CHEBFUN F.
%
%   SINH(F, PREF) does the same but uses the CHEBPREF object PREF when computing
%   the composition.
%
% See also ASINH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Call the compose method:
F = compose(F, @sing, pref);

end

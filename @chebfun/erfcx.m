function g = erfcx(f, pref)
%ERFCX   Scaled complementary error function of a CHEBFUN.
%   Y = ERFCX(X) is the scaled complementary error function of the CHEBFUN X.
%   X must be real.  The scaled complementary error function is defined as:
%       ERFCX(X) = EXP(X.^2) * ERFC(X)
%   which is approximately (1/sqrt(pi)) * 1./X for large X.
%
% See also ERF, ERFC, ERFINV, ERFCINV.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Input must be real:
if ( ~isreal(f) )
    error('CHEBFUN:erfcx:notreal', 'Input must be real.');
end

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Call the compose method:
g = compose(f, @erfcx, pref);

end

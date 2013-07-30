function g = erfc(f, pref)
%ERFC   Complementary error function of a chebfun.
%   Y = ERFC(X) is the complementary error function for the chebfun X. X must be
%   real. The complementary error function is defined as:
%       ERFC(X)(s) = 2/sqrt(pi) * integral from X(s) to inf of exp(-t^2) dt.
%               = 1 - ERF(X)(s).
%
% See also ERF, ERFCX, ERFINV, ERFCINV.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Input must be real:
if ( ~isreal(f) )
    error('CHEBFUN:erfc:notreal', 'Input must be real.');
end

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @erfc, pref);

end
function g = erf(f, pref)
%ERF   Error function of a chebfun.
%   Y = ERF(X) is the error function for each element of X. X must be real.
%
%   The error function is defined as:
%       erf(X)(s) = 2/sqrt(pi) * integral from 0 to X(s) of exp(-t^2) dt.
%
% See also ERFC, ERFCX, ERFINV, ERFCINV.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Input must be real:
if ( ~isreal(f) )
    error('CHEBFUN:erf:notreal', 'Input must be real.');
end

% Obtain preferences:
if ( nargin == 1 )
    pref = chebfun.pref();
end

% Call the compose method:
g = compose(f, @erf, pref);

end
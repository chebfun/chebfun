function F = erf(F, pref)
%ERF   Error function of a CHEBFUN.
%   ERF(X) is the error function of the CHEBFUN X. X must be real.
%
%   The error function is defined as:
%       erf(X)(s) = 2/sqrt(pi) * integral from 0 to X(s) of exp(-t^2) dt.
%
% See also ERFC, ERFCX, ERFINV, ERFCINV.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Input must be real:
if ( ~isreal(F) )
    error('CHEBFUN:erf:notreal', 'Input must be real.');
end

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Loop over the columns of F:
for k = 1:numel(F)
    % Call the compose method:
    F(k) = compose(F(k), @erf, pref);
end

end

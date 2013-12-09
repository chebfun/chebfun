function F = erfc(F, pref)
%ERFC   Complementary error function of a CHEBFUN.
%   ERFC(X) is the complementary error function for the CHEBFUN X. X must be
%   real. The complementary error function is defined as:
%       ERFC(X)(s) = 2/sqrt(pi) * integral from X(s) to inf of exp(-t^2) dt.
%                  = 1 - ERF(X)(s).
%
% See also ERF, ERFCX, ERFINV, ERFCINV.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Input must be real:
if ( ~isreal(F) )
    error('CHEBFUN:erfc:notreal', 'Input must be real.');
end

% Obtain preferences:
if ( nargin == 1 )
    pref = chebpref();
end

% Loop over the columns of F:
for k = 1:numel(F)
    % Call the compose method:
    F(k) = compose(F(k), @erfc, pref);
end

end

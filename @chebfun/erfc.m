function F = erfc(F, varargin)
%ERFC   Complementary error function of a CHEBFUN.
%   ERFC(X) is the complementary error function of the real-valued CHEBFUN X.
%   The complementary error function is defined as:
%       ERFC(X)(s) = 2/sqrt(pi) * integral from X(s) to inf of exp(-t^2) dt.
%                  = 1 - ERF(X)(s).
%
% See also ERF, ERFCX, ERFINV, ERFCINV.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Input must be real:
if ( ~isreal(F) )
    error('CHEBFUN:CHEBFUN:erfc:notreal', 'Input must be real.');
end

% Call the compose method:
F = compose(F, @erfc, varargin{:});

end

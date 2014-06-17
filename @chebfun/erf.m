function F = erf(F, varargin)
%ERF   Error function of a CHEBFUN.
%   ERF(X) is the error function of the real-valued CHEBFUN X.
%
%   The error function is defined as:
%       erf(X)(s) = 2/sqrt(pi) * integral from 0 to X(s) of exp(-t^2) dt.
%
% See also ERFC, ERFCX, ERFINV, ERFCINV.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Input must be real:
if ( ~isreal(F) )
    error('CHEBFUN:CHEBFUN:erf:notreal', 'Input must be real.');
end

% Call the compose method:
F = compose(F, @erf, varargin{:});

end

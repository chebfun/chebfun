function F = erfcinv(F, varargin)
%ERFCINV   Inverse complementary error function of a CHEBFUN.
%   X = ERFCINV(Y) is the inverse of the complementary error function for the
%   CHEBFUN Y. It satisfies Y = ERFC(X) for 2 >= Y >= 0 and -Inf <= X <= Inf.
%
% See also ERF, ERFC, ERFCX, ERFINV.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Input must be real:
if ( ~isreal(F) )
    error('CHEBFUN:CHEBFUN:erfcinv:notreal', 'Input must be real.');
end

% Call the compose method:
F = compose(F, @erfcinv, varargin{:});

end

function F = erfcx(F, varargin)
%ERFCX   Scaled complementary error function of a CHEBFUN.
%   ERFCX(X) is the scaled complementary error function of the real-valued
%   CHEBFUN X. The scaled complementary error function is defined as:
%       ERFCX(X) = EXP(X.^2) * ERFC(X)
%   which is approximately (1/sqrt(pi)) * 1./X for large X.
%
% See also ERF, ERFC, ERFINV, ERFCINV.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Call the compose method:
F = compose(F, @erfcx, varargin{:});

end

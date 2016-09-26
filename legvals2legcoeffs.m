function legcoeffs = legvals2legcoeffs(legvals)
%LEGVALS2LEGCOEFFS  Convert Legendre values to Legendre coefficients.
% 	LEGCOEFFS = LEGVALS2LEGCOEFFS(LEGVALS), converts the column vector
%   LEGVALS representing values on a Legendre grid (i.e, F(LEGPTS)) to a
%   vector LEGCOEFFS representing the Chebyshev coefficients of the series
%       F(X) = C_LEG(1)*P0(X) + ... + C_LEG(N)*P{N-1}(X).
% 
% See also CHEBFUN.IDLT, LEGPTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% This command is a wrapper for chebfun/idlt.
legcoeffs = chebfun.idlt(legvals);

end

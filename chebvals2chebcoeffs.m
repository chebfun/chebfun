function chebcoeffs = chebvals2chebcoeffs( chebvals )
%CHEBVALS2CHEBCOEFFS  Convert Chebyshev values to coefficients.
% 	CHEBCOEFFS = CHEBVALS2CHEBCOEFFS(CHEBVALS), converts the column vector
%   CHEBVALS representing values on a Chebyshev grid (i.e, F(CHEBPTS)) to a
%   vector CHEBCOEFFS representing the Chebyshev coefficients of the series
%       F(X) = C_CHEB(1)*T0(X) + ... + C_CHEB(N)*T{N-1}(X).
% 
% See also CHEBTECH2.VALS2COEFFS, CHEBPTS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% This command is a wrapper for chebtech2/vals2coeffs.
chebcoeffs = chebtech2.vals2coeffs( chebvals ); 

end
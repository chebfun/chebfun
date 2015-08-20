function chebvals = chebcoeffs2chebvals( chebcoeffs )
%CHEBCOEFFS2CHEBVALS  Convert Chebyshev coefficients to values.
%   CHEBVALS = CHEBCOEFFS2CHEBVALS(CHEBCOEFFS), converts the columns vector
%   CHEBCOEFFS representing Chebyshev coefficients in the Chebyshev series
%       F(X) = C_CHEB(1)*T0(X) + ... + C_CHEB(N)*T{N-1}(X), 
%   to a vector CHEBVALS representing values of the expansion at CHEBPTS, i.e., 
%       CHEBVALS = F(chebpts(N));
% 
% See also CHEBTECH2.COEFFS2VALS, CHEBPTS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% This command is a wrapper for chebtech2/coeffs2vals:
chebvals = chebtech2.coeffs2vals( chebcoeffs ); 

end
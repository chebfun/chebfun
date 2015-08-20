function legvals = chebvals2legvals( chebvals )
% CHEBVALS2LEGVALS      Convert Chebyshev values to Legendre values
% 
%  LEGVALS = CHEBVALS2LEGVALS( CHEBVALS ), converts a vector of
%  CHEBVALS representing values of a Chebyshev expansion at CHEBPTS to a
%  vector LEGVALS representing the expansion evaluated at LEGPTS.  
% 
% See also vals2coeffs, chebpts.

chebcoeffs = chebtech2.vals2coeffs( chebvals ); 
legvals = chebfun.ndct( chebcoeffs ); 

end
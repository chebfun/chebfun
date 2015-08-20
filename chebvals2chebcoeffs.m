function chebcoeffs = chebvals2chebcoeffs( chebvals )
% CHEBVALS2CHEBCOEFFS       Convert Chebyshev values to coefficients
% 
%  CHEBCOEFFS = CHEBVALS2CHEBCOEFFS( CHEBVALS ), converts the vector of
%  CHEBVALS representing values of a Chebyshev expansion at CHEBPTS to a
%  vector CHEBCOEFFS representing the Chebyshev coefficients of the 
%  expansion.  
%
%  This command is a wrapper for chebtech2/vals2coeffs.
% 
% See also vals2coeffs, chebpts.

chebcoeffs = chebtech2.vals2coeffs( chebvals ); 

end
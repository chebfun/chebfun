function legcoeffs = legvals2legcoeffs( legvals )
% LEGVALS2LEGCOEFFS       Convert Legendre values to Legendre coeffs
% 
%  LEGCOEFFS = LEGVALS2CHEBCOEFFS( LEGVALS ), converts the vector of
%  LEGVALS representing values of a Legendre expansion at LEGPTS to a
%  vector LEGCOEFFS representing the Legendre coefficients of the 
%  expansion.  
%
%  This command is a wrapper for chebfun/idlt.
% 
% See also idlt.

legcoeffs = chebfun.idlt( legvals ); 

end
function legvals = legcoeffs2legvals( legcoeffs )
% LEGCOEFFS2LEGVALS       Convert Legendre coeffs to Legendre values
% 
%  LEGVALS = LEGCOEFFS2LEGVALS( LEGCOEFFS ), converts the vector of
%  LEGCOEFFS representing Legendre coefficients in an expansion to a vector 
%  vector LEGVALS representing values of the expansion at LEGPTS.  
%
%  This command is a wrapper for chebfun/dlt.
% 
% See also dlt.

legvals = chebfun.dlt( legcoeffs ); 

end
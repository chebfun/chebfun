function chebvals = legvals2chebvals( legvals )
% LEGVALS2CHEBVALS       Convert Legendre values to Chebyshev values
% 
%  CHEBVALS = LEGVALS2CHEBVALS( LEGVALS ), converts the vector of
%  LEGVALS representing values of a polynomial at LEGPTS to a
%  vector CHEBVALS representing the same polynomial evaluated at CHEBPTS.  
% 
% See also chebvals2legvals

legcoeffs = chebfun.idlt( legvals ); 
chebvals = legcoeffs2chebvals( legcoeffs ); 

end
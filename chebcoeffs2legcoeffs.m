function legcoeffs = chebcoeffs2legcoeffs( chebcoeffs )
% CHEBCOEFFS2LEGCOEFFS       Convert Chebyshev coeffs to Legendre coeffs
% 
%  LEGCOEFFS = CHEBCOEFFS2LEGCOEFFS( CHEBCOEFFS ), converts a vector of
%  CHEBCOEFFS representing Chebyshev coefficients of an expansion to 
%  a vector LEGCOEFFS representing Legendre coefficients. That is, 
% 
%     sum_k  CHEBCOEFFS( k ) T_k(x)    =  sum_k LEGCOEFFS( k )  P_k(x)
%
%  This command is a wrapper for chebfun/cheb2leg.
% 
% See also cheb2leg.

legcoeffs = chebfun.cheb2leg( chebcoeffs ); 

end
function chebcoeffs = legcoeffs2chebcoeffs( legcoeffs )
% LEGCOEFFS2CHEBCOEFFS      Convert Legendre coeffs to Chebyshev coeffs
% 
%  CHEBCOEFFS = LEGCOEFFS2CHEBCOEFFS( LEGCOEFFS ), converts a vector of
%  LEGCOEFFS representing Legendre coefficients of an expansion to 
%  a vector CHEBCOEFFS representing Chebyshev coefficients. That is,   
% 
%      sum_k  LEGCOEFFS( k ) P_k(x)    =  sum_k CHEBCOEFFS( k )  T_k(x)
% 
% See also vals2coeffs, chebpts.

chebcoeffs = chebfun.leg2cheb( legcoeffs ); 

end
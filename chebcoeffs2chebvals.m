function chebvals = chebcoeffs2chebvals( chebcoeffs )
% CHEBCOEFFS2CHEBVALS       Convert Chebyshev coefficients to values 
% 
%  CHEBVALS = CHEBCOEFFS2CHEBVALS( CHEBCOEFFS ), converts the vector of
%  CHEBCOEFFS representing Chebyshev coefficients in an expansion to 
%  a vector CHEBVALS representing values of the expansion at CHEBPTS.
%
%  This command is a wrapper for chebtech2/coeffs2vals.
% 
% See also coeffs2vals, chebpts.

chebvals = chebtech2.coeffs2vals( chebcoeffs ); 

end
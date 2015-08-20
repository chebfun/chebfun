function legvals = chebvals2legvals( chebvals )

chebcoeffs = chebtech2.vals2coeffs( chebvals ); 
legvals = chebfun.ndct( chebcoeffs ); 

end
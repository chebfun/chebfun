function legcoeffs = legvals2legcoeffs( legvals )

legcoeffs = chebfun.idlt( legvals ); 

end
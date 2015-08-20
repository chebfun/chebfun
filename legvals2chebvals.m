function chebvals = legvals2chebvals( legvals )

legcoeffs = chebfun.idlt( legvals ); 
chebvals = legcoeffs2chebvals( legcoeffs ); 

end
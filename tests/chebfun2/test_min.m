function pass = test_min( pref ) 
% Test the chebfun2/min command. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1000*pref.eps; 
tol = sqrt(tol); 

f = chebfun2(@(x,y) cos(x.*y)); 
g = chebfun(@(x) cos(x)); 

h1 = min( f ); 
h2 = min( f, []); 
h3 = min( f, [], 2 );
h4 = min( f, [], 3 );  

pass(1) = ( norm( h1 - g.' ) < tol );
pass(2) = ( norm( h2 - g.' ) < tol );
pass(3) = ( norm( h3 - g ) < tol );
pass(4) = ( norm( h4 - f ) < tol );


f = chebfun2(@(x,y) x, [-2 3 -4 10]); 

h1 = min( f ); 
h2 = min( f, []); 
h3 = min( f, [], 2 );
h4 = min( f, [], 3 );  

pass(5) = ( norm( h1 - chebfun(@(x) x,[-2,3]).' ) < tol );
pass(6) = ( norm( h2 - chebfun(@(x) x,[-2,3]).' ) < tol );
pass(7) = ( norm( h3 - chebfun(@(x) -2+0*x,[-4,10]) ) < 10*tol );
pass(8) = ( norm( h4 - f ) < tol );
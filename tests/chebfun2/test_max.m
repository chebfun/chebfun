function pass = test_max( pref ) 
% Test the chebfun2/max command. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 100*pref.eps; 

f = chebfun2(@(x,y) cos(x.*y)); 
g = chebfun(@(x) 1 + 0*x); 

h1 = max( f ); 
h2 = max( f, []); 
h3 = max( f, [], 2 );
h4 = max( f, [], 3 );  

pass(1) = ( norm( h1 - g.' ) < tol );
pass(2) = ( norm( h2 - g.' ) < tol );
pass(3) = ( norm( h3 - g ) < tol );
pass(4) = ( norm( h4 - f ) < tol );


f = chebfun2(@(x,y) x, [-2 3 -4 10]); 

h1 = max( f ); 
h2 = max( f, []); 
h3 = max( f, [], 2 );
h4 = max( f, [], 3 );  

pass(5) = ( norm( h1 - chebfun(@(x) x,[-2,3]).' ) < tol );
pass(6) = ( norm( h2 - chebfun(@(x) x,[-2,3]).' ) < tol );
pass(7) = ( norm( h3 - chebfun(@(x) 3+0*x,[-4,10]) ) < 1e12*tol );
pass(8) = ( norm( h4 - f ) < tol );
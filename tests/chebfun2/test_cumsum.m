function pass = test_cumsum( pref )

if ( nargin == 0 )
    pref = chebpref; 
end

tol = 100*pref.cheb2Prefs.eps; 
j = 1; 

x = chebfun2(@(x,y) x, [-1 1 -1 1]); 
y = chebfun2(@(x,y) y, [-1 1 -1 1]); 

pass(j) = ( norm( cumsum(x) - x.*(y+1) ) < tol ); j = j + 1;
pass(j) = ( norm( cumsum(x, 2) - .5*(x.^2-1) ) < tol ); j = j + 1;
pass(j) = ( norm( cumsum(y, 2) - y.*(x+1) ) < tol ); j = j + 1;
pass(j) = ( norm( cumsum(y) - .5*(y.^2-1) ) < tol ); j = j + 1;

x = chebfun2(@(x,y) x, [-1.1 2 -.2 3]); 
y = chebfun2(@(x,y) y, [-1.1 2 -.2 3]); 

pass(j) = ( norm( cumsum(x) - x.*(y+.2) ) < tol ); j = j + 1;
pass(j) = ( norm( cumsum(x, 2) - .5*(x.^2-1.1^2) ) < tol ); j = j + 1;
pass(j) = ( norm( cumsum(y, 2) - y.*(x+1.1) ) < tol ); j = j + 1;
pass(j) = ( norm( cumsum(y) - .5*(y.^2-.2^2) ) < tol ); j = j + 1;


end
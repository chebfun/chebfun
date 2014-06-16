function pass = test_cumsum( pref )

if ( nargin == 0 )
    pref = chebfunpref; 
end

tol = 100*pref.eps; 
j = 1; 

% Check cumsum on square domain
x = chebfun2(@(x,y) x, [-1 1 -1 1]); 
y = chebfun2(@(x,y) y, [-1 1 -1 1]); 

pass(j) = ( norm( cumsum(x) - x.*(y+1) ) < tol ); j = j + 1;
pass(j) = ( norm( cumsum(x, 2) - .5*(x.^2-1) ) < tol ); j = j + 1;
pass(j) = ( norm( cumsum(y, 2) - y.*(x+1) ) < tol ); j = j + 1;
pass(j) = ( norm( cumsum(y) - .5*(y.^2-1) ) < tol ); j = j + 1;

% Check cumsum on rectangular domain
x = chebfun2(@(x,y) x, [-1.1 2 -.2 3]); 
y = chebfun2(@(x,y) y, [-1.1 2 -.2 3]); 

pass(j) = ( norm( cumsum(x) - x.*(y+.2) ) < tol ); j = j + 1;
pass(j) = ( norm( cumsum(x, 2) - .5*(x.^2-1.1^2) ) < tol ); j = j + 1;
pass(j) = ( norm( cumsum(y, 2) - y.*(x+1.1) ) < tol ); j = j + 1;
pass(j) = ( norm( cumsum(y) - .5*(y.^2-.2^2) ) < tol ); j = j + 1;

% Check that a double cumsum is a cumsum2 
f = sin((x-.1).*(y+.1)); 
pass(j) = ( norm( cumsum(cumsum(f),2) - cumsum2(f) ) < tol ); j = j + 1;
pass(j) = ( norm( cumsum(cumsum(f,2)) - cumsum2(f) ) < tol ); j = j + 1;
end
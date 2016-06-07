function pass = test_lu( pref ) 
% Test for LU decomposition of a chebfun2. 

if ( nargin == 0 ) 
    pref = chebfunpref; 
end 

tol = 100*pref.cheb2Prefs.chebfun2eps;

% Decomposition on [-1,1,-1,1]: 
x = chebpts(100);
f = chebfun2( @(x,y) cos( x.*y ) );
[L, U] = lu( f );
g = L * U; 
pivPos = f.pivotLocations; 

% Accurate decomposition: 
pass(1) = norm( fevalm(f, x, x) - fevalm(g, x, x) ) < tol; 
pass(2) = norm( diag( L( pivPos( :, 2 ) , :) ) - ones( length(f), 1) ) < tol;
% pass(j) = norm( triu( L( pivPos( :, 2 ) , :) , 1) ) < sqrt(tol); j = j + 1; 
pass(3) = norm( tril( U( :, pivPos( :, 1 ) ) , -1) ) < sqrt(tol);



% Try the same thing on different domain: 
x = chebpts(100, [-2.1, 4.3]);
y = chebpts(100, [-1, 2.7]);

f = chebfun2( @(x,y) cos( x.*y ), [-2.1, 4.3, -1, 2.7]);
[L, U] = lu( f );
g = L * U; 
pivPos = f.pivotLocations; 

% Accurate decomposition: 
pass(4) = norm( fevalm(f, x, y) - fevalm(g, x, y) ) < 2*tol; 
pass(5) = norm( diag( L( pivPos( :, 2 ) , :) ) - ones( length(f), 1) ) < tol; 
% pass(j) = norm( triu( L( pivPos( :, 2 ) , :) , 1) ) < sqrt(tol); j = j + 1; 
pass(6) = norm( tril( U( :, pivPos( :, 1 ) ) , -1) ) < sqrt(tol); 


end
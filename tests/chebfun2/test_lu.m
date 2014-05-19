function pass = test_lu( pref ) 
% Test for LU decomposition of a chebfun2. 

if ( nargin == 0 ) 
    pref = chebfunpref; 
end 

tol = 100*pref.cheb2Prefs.eps; 
j = 1; 

% Decomposition on [-1,1,-1,1]: 
x = chebpts(100);
f = chebfun2( @(x,y) cos( x.*y ) );
[L, U] = lu( f );
g = L * U; 
pivPos = f.pivotLocations; 

% Accurate decomposition: 
pass(j) = norm( fevalm(f, x, x) - fevalm(g, x, x) ) < tol; j = j + 1; 
pass(j) = norm( diag( L( pivPos( :, 2 ) , :) ) - ones( length(f), 1) ) < tol; j = j + 1; 
% pass(j) = norm( triu( L( pivPos( :, 2 ) , :) , 1) ) < sqrt(tol); j = j + 1; 
pass(j) = norm( tril( U( :, pivPos( :, 1 ) ) , -1) ) < sqrt(tol); j = j + 1; 



% Try the same thing on different domain: 
x = chebpts(100, [-2.1, 4.3]);
y = chebpts(100, [-1, 2.7]);

f = chebfun2( @(x,y) cos( x.*y ), [-2.1, 4.3, -1, 2.7]);
[L, U] = lu( f );
g = L * U; 
pivPos = f.pivotLocations; 

% Accurate decomposition: 
pass(j) = norm( fevalm(f, x, y) - fevalm(g, x, y) ) < tol; j = j + 1; 
pass(j) = norm( diag( L( pivPos( :, 2 ) , :) ) - ones( length(f), 1) ) < tol; j = j + 1; 
% pass(j) = norm( triu( L( pivPos( :, 2 ) , :) , 1) ) < sqrt(tol); j = j + 1; 
pass(j) = norm( tril( U( :, pivPos( :, 1 ) ) , -1) ) < sqrt(tol); j = j + 1; 


end
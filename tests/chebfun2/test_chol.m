function pass = test_chol( pref ) 
% Test for Cholesky decomposition of a chebfun2. 

if ( nargin == 0 ) 
    pref = chebfunpref; 
end 

tol = 100*pref.eps; 
j = 1; 

% Decomposition on [-1,1,-1,1]: 
x = chebpts(100);
f = chebfun2( @(x,y) 1 ./ ( 10 + x.^2 + y.^2 ) );
R = chol( f );
g = R' * R; 
pivPos = f.pivotLocations; 

% Accurate decomposition: 
pass(j) = norm( pivPos(:, 1) - pivPos(:, 2) ) < tol; j = j + 1; 
pass(j) = norm( fevalm(f, x, x) - fevalm(g, x, x) ) < tol; j = j + 1; 
pass(j) = norm( R( 1, pivPos( 1, 2 ) ) - sqrt(feval(f, pivPos(1,1), pivPos(1,2))) ) < tol; j = j + 1; 
pass(j) = norm( tril( R( :, pivPos( :, 2 ) ) , -1) ) < tol; j = j + 1; 

%%

% Try the same thing on different domain: 
x = chebpts(100, [-2.1 4.3]);

f = chebfun2( @(x,y) 1 ./ ( 10 + x.^2 + y.^2 ) , [-2.1 4.3 -2.1 4.3]);
R = chol( f );
g = R' * R; 
pivPos = f.pivotLocations; 

% Accurate decomposition: 
pass(j) = norm( pivPos(:, 1) - pivPos(:, 2) ) < tol; j = j + 1; 
pass(j) = norm( fevalm(f, x, x) - fevalm(g, x, x) ) < tol; j = j + 1; 
pass(j) = norm( R( 1, pivPos( 1, 2 ) ) - sqrt(feval(f, pivPos(1,1), pivPos(1,2))) ) < tol; j = j + 1; 
pass(j) = norm( tril( R( :, pivPos( :, 2 ) ) , -1) ) < tol; j = j + 1; 


% KX's bug for chebfun2 objects that are constructed from data: 
x = ones(1,2); 
f = chebfun2(x,'coeffs'); 
g = chol(f);
% g should be linear and a row: 
pass(j) = norm( g - chebfun([0;g(1)]).' ) < tol; j = j + 1; 
end
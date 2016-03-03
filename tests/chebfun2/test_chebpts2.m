function pass = test_chebpts2( pref ) 
% Test CHEBPTS2 

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 10*pref.cheb2Prefs.chebfun2eps;

% One arguments: 
n = 10; 
d = [-1, 1];
[xx1, yy1] = chebfun2.chebpts2(n); 
x = chebpts( n, d ); 
y = chebpts( n, d ); 
[xx2, yy2] = meshgrid( x, y ); 
pass(1) = ( norm( xx1 - xx2 ) < tol ); 
pass(2) = ( norm( yy1 - yy2 ) < tol ); 

% Two arguments:
n = 10; 
m = 7; 
d = [-1, 1];
[xx1, yy1] = chebfun2.chebpts2(n, m); 
x = chebpts( n, d); 
y = chebpts( m, d ); 
[xx2, yy2] = meshgrid( x, y ); 
pass(3) = ( norm( xx1 - xx2 ) < tol ); 
pass(4) = ( norm( yy1 - yy2 ) < tol );

% Three arguments:
n = 10; 
m = 7; 
[xx1, yy1] = chebfun2.chebpts2(n, m, [-2 1 -1 2]); 
x = chebpts( n, [-2 1]); 
y = chebpts( m, [-1 2] ); 
[xx2, yy2] = meshgrid( x, y ); 
pass(5) = ( norm( xx1 - xx2 ) < 10*tol ); 
pass(6) = ( norm( yy1 - yy2 ) < tol );

end
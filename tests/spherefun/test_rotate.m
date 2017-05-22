function pass = test_rotate( pref ) 
% Check that roots works for a spherefun.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e4*pref.cheb2Prefs.chebfun2eps;

% Rotate a random function: 
f = randnfunsphere( 0.05 );
g = rotate(f, 0.1, 0.1, 0.1);
h = rotate(g, -0.1, -0.1, -0.1);
pass(1) = ( norm( f - h ) < tol );

end
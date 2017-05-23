function pass = test_rotate( pref ) 
% Check that roots works for a spherefun.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e3*pref.cheb2Prefs.chebfun2eps;

% Rotate a random function: 
angls = [0.1 0.1 0.1];
f = randnfunsphere( 0.05 );
g = rotate(f, angls(1), angls(2), angls(3));
% Undo the rotation
h = rotate(g, -angls(3), -angls(2), -angls(1));
pass(1) = ( norm( f - h ) < 10*tol );

% Random choice of rotation angles
angls = [0.704147450733711 1.158707918472494 0.884902035498222];
g = rotate(f, angls(1), angls(2), angls(3));
% Undo the rotation
h = rotate(g, -angls(3), -angls(2), -angls(1));
pass(2) = ( norm( f - h ) < 10*tol );

% Check the that the two methods for doing the rotation produce the same
% result
f = spherefun(@(lam,th) sin(cos(th) + cos(lam-0.2).*sin(th) + sin(lam+0.4).*sin(th)).^8);
g = rotate(f, angls(1), angls(2), angls(3), 'feval');
h = rotate(f, angls(1), angls(2), angls(3), 'nufft');
pass(3) = ( norm( g - h ) < 10*tol );

% Check that the integral of a rotated spherefun is preserved.
g = rotate(f, angls(1), angls(2), angls(3));
pass(4) = ( abs( sum2(g) - sum2(f) ) < tol );

% Rotate spherical harmonics with various symmetries:
f = spherefun.sphharm(21, 0);  % Only a function of z
g = rotate(f, sqrt(2), 0, 0);
pass(5) = ( norm( f - g ) < tol );

f = spherefun.sphharm(10, 10); 
g = rotate(f, pi/5, 0, 0);
pass(6) = ( norm( f - g ) < tol );

% Function with octahedral symmetry
f = spherefun.sphharm(4, 0) + sqrt(5/7)*spherefun.sphharm(4, 4);
g = rotate(f, pi/2, pi/2, 0);
pass(7) = ( norm( f - g ) < tol );
g = rotate(g, pi/2, pi/2, 0);
pass(8) = ( norm( f - g ) < tol );

end
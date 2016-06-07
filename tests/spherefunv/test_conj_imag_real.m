function pass = test_conj_imag_real( pref ) 
% Test CONJ, IMAG, REAL

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 100*pref.cheb2Prefs.chebfun2eps;

% Test function
f = spherefun(@(x,y,z) cos((x+.1).*y.*z));
% Spherefunv
u = grad(f);

% Get some random points on the sphere for testing.
rng(7); lam0 = rand; th0 = rand;

% u is real so conjugate is real and same as u
v = conj( u ); 
pass(1) = ( norm( u(lam0,th0) - v(lam0,th0) ) < tol ); 

% u is real so real(u) is the same as u
v = real( u );
pass(2) = ( norm( u(lam0,th0) - v(lam0,th0) ) < tol ); 

% u is real so imag(u) is zero
v = imag( u );
pass(3) = ( norm( v ) < tol ); 

% Developer note: If complex-valued spherefun and spherefunv objects are
% allowed then more tests should be added here.

end
function pass = test_conj_imag_real( pref ) 
% Test CONJ, IMAG, REAL

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 100*pref.cheb2Prefs.chebfun2eps;

% Test function
f = diskfun(@(x,y) cos((x+.1).*y));
% diskfunv
u = grad(f);
% Get some random points on the disk for testing.
rng(7); r0 = rand; th0 = pi*(2*rand-1);

% u is real so conjugate is real and same as u
v = conj( u ); 
pass(1) = ( norm( u(th0,r0, 'polar') - v(th0,r0, 'polar') ) < tol ); 

% u is real so real(u) is the same as u
v = real( u );
pass(2) = ( norm( u(th0,r0, 'polar') - v(th0,r0, 'polar') ) < tol ); 

% u is real so imag(u) is zero
v = imag( u );
pass(3) = ( norm( v ) < tol ); 

% Developer note: If complex-valued diskfun and diskfunv objects are
% allowed then more tests should be added here.

end
function pass = test_conj(pref)
% Test CONJ

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 100*pref.cheb3Prefs.chebfun3eps;

ff = @(x,y,z) cos(x.*y.*z);

f = chebfun3v(ff, ff); 
g = conj(f); 
pass(1) = norm(f - g) < tol;

f = chebfun3v(ff,ff); 
g = conj( 1i*f ); 
pass(2) = norm(1i*f + g) < tol;

f1 = chebfun3v(ff, ff);
f2 = chebfun3v(@(x,y,z) sin(x + y.^2), @(x,y,z) sin(x + y.^2));
g = conj(f1 + 1i*f2); 
pass(3) = norm(f1 - 1i*f2 - g) < tol;

end
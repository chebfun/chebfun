function pass = test_real(pref)
% Test REAL

if ( nargin == 0 ) 
    pref = chebfunpref;
end
tol = 100*pref.cheb3Prefs.chebfun3eps;

f = chebfun3v(@(x,y,z) cos(x.*y.*z), @(x,y,z) cos(x.*y.*z));
g = real(f); 
pass(1) = norm(g - f) < tol;

f = chebfun3v(@(x,y,z) cos(x.*y.*z), @(x,y,z) cos(x.*y.*z));
g = real(1i*f);
pass(2) = norm(g) < tol;

f1 = chebfun3v(@(x,y,z) cos(x.*y.*z), @(x,y,z) cos(x.*y.*z));
f2 = chebfun3v(@(x,y,z) sin(x + y.^2+z.^3), @(x,y,z) sin(x + y.^2+z.^3));
g = real(f1 + 1i*f2);
pass(3) = norm(g - f1) < tol;

end
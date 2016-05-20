function pass = test_integral(pref)
% Test path integral for 3D functions.

if ( nargin == 0)
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) 2*x.*y+3*z); 
curve = chebfun(@(t) [2*t 5*t 4*t], [0, 1]); 
exact = 38*sqrt(5);
pass(1) = abs(integral(f, curve) - exact) < tol;

f = chebfun3(@(x,y,z) x.^2.*z); 
curve = chebfun(@(t) [-t 6+2*t 2+5*t], [0 1]);
exact = 23*sqrt(30)/12;
pass(2) = abs(integral(f, curve) - exact ) < tol;

end
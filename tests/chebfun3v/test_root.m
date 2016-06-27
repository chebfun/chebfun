function pass = test_root(pref)
% Check finding one common root of a CHEBFUN3V object.

if ( nargin < 1 ) 
    pref = chebfunpref;
end 
tol = 1e3 * pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) y-x.^2);
g = chebfun3(@(x,y,z) z-x.^3);
h = chebfun3(@(x,y,z) cos(exp(x.*sin(-2+y+z))));
F = [f; g; h];
r = root(F);

pass(1) = numel(r) == 3;  

pass(2) = all(F(r(1), r(2), r(3)) < [tol; tol; tol]);
end
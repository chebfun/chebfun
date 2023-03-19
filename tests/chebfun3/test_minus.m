function pass = test_minus(pref)
% Test for chebfun3/minus.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.cheb3Prefs.chebfun3eps;

ff = @(x,y,z) cos(x.*y.*z);
gg = @(x,y,z) x + y + z + x.*y.*z;

f = chebfun3(ff);
g = chebfun3(gg);
FminusG = chebfun3(@(x,y,z) ff(x,y,z) - gg(x, y, z));
pass(1) = norm((f-g) - FminusG) < tol;

dom = [-1 pi 0 2*pi -pi pi];
f = chebfun3(ff, dom);
g = chebfun3(gg, dom);
FminusG = chebfun3(@(x,y,z) ff(x,y,z) - gg(x, y, z), dom);
pass(2) = norm((f-g) - FminusG) < tol;

[x, ~, ~] = cheb.xyz;
f = x;
%g = 0.9999999999*x;
g = 0.9999999*x;
% f - g should NOT flush to zero. This makes sure that the tolerance set in
% chebfun3/plus does not hugely overestimate.
pass(3) = norm(f-g) > 0;

end
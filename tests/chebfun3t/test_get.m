function pass = test_get(pref)
% Test CHEBFUN3T/GET.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.cheb3Prefs.chebfun3eps;

f = chebfun3t(@(x,y,z) cos(x.*y.*z)); 
pass(1) = norm([-1 1 -1 1 -1 1] - f.domain) < tol; 

pass(2) = f.vscale > 0;

pass(3) = norm(f.coeffs(:)) > 0;
end
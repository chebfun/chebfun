function pass = test_get(pref)
% Test CHEBFUN3T/GET.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.cheb3Prefs.chebfun3eps;
j = 1; 

f = chebfun3t(@(x,y,z) cos(x.*y.*z)); 
pass(j) = norm([-1 1 -1 1 -1 1] - f.domain) < tol; 
j = j + 1; 

pass(j) = f.vscale > 0; 
j = j + 1; 

pass(j) = norm(f.coeffs(:)) > 0; 
j = j + 1; 
end
function pass = test_mean2(pref)
% Test mean2 command of Chebfun3.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) x.^2.*y.^4.*exp(z));
g = chebfun(@(z) exp(z)/15);
pass(1) = norm(mean2(f) - g) < tol; 

end
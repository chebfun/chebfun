function pass = test_mean(pref)
% Test mean command of Chebfun3.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) x.^2.*y.^4.*exp(z));
g = chebfun2(@(y,z) y.^4.*exp(z)/3);
pass(1) = norm(mean(f) - g) < tol;

end
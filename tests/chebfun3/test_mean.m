function pass = test_mean(pref)
% Test mean3, mean2 and mean commands of Chebfun3.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) x.^2.*y.^4.*exp(z));
exact = (exp(1)-exp(-1))/30;
pass(1) = abs(mean3(f) - exact) < tol; 

g = chebfun(@(z) exp(z)/15);
pass(2) = norm(mean2(f) - g) < tol; 

g = chebfun2(@(y,z) y.^4.*exp(z)/3);
pass(3) = norm(mean(f) - g) < tol;

end
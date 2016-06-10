function pass = test_mean( pref )
% Test mean3, mean2 and mean commands of Chebfun3.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 

tol = pref.cheb3Prefs.chebfun3eps;
j = 1; 

f = chebfun3(@(x,y,z) x.^2.*y.^4.*exp(z));
exact = (exp(1)-exp(-1))/30;
pass(j) = abs(mean3(f) - exact) < tol; 
j = j + 1; 

g = chebfun(@(z) exp(z)/15);
pass(j) = norm(mean2(f) - g) < tol; 
j = j + 1; 

g = chebfun2(@(y,z) y.^4.*exp(z)/3);
pass(j) = norm(mean(f) - g) < tol; 
j = j + 1;

end
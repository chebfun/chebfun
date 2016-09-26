function pass = test_mean3(pref)
% Test mean3 command of Chebfun3.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 100*pref.cheb3Prefs.chebfun3eps;

% First test:
f = chebfun3(@(x,y,z) x.^2.*y.^4.*exp(z));
exact = (exp(1)-exp(-1))/30;
pass(1) = abs(mean3(f) - exact) < tol; 

% A second test:
ff = @(x,y,z) sin(pi*x).^2+ sin(pi*(y+z)).^2;
f = chebfun3(ff);
exact = 1;
pass(2) = abs(mean3(f) - exact) < tol;

end
function pass = test_imag(pref) 
% Test IMAG

if ( nargin == 0 ) 
    pref = chebfunpref; 
end
tol = 10*pref.cheb3Prefs.chebfun3eps;

f = chebfun3v(@(x,y,z) cos(x.*y.*z), @(x,y,z) sin(x.*y.*z), ...
    @(x,y,z) exp(x.*y.*z));
g = imag(f); 
pass(1) = norm(g) < tol;

f = chebfun3v(@(x,y,z) cos(x.*y.*z), @(x,y,z) sin(x.*y.*z), ...
    @(x,y,z) exp(x.*y.*z));
g = imag(1i*f);
pass(2) = norm(f - g) < 100*tol;

f1 = chebfun3v(@(x,y,z) cos(x.*y.*z), @(x,y,z) sin(x.*y.*z), ...
    @(x,y,z) exp(x.*y.*z));
f2 = chebfun3v(@(x,y,z) sin(x + y.^2 + z), @(x,y,z) sin(x + y.^2 + z),...
    @(x,y,z) exp(x.*y.*z));
g = imag(f1 + 1i*f2);
pass(3) = norm(f2 - g) < 100*tol;

end
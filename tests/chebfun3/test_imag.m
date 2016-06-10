function pass = test_imag( pref ) 
% Test IMAG
if ( nargin == 0 ) 
    pref = chebfunpref; 
end
tol = 100*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
g = chebfun3(@(x,y,z) sin(x+y.^2+z.^3));

% Simple consistency check: 
x = linspace(-1, 1, 3);
[xx, yy, zz] = ndgrid(x);
h = f + 1i*g;
fVals = feval(h, xx, yy, zz);
gVals = feval(g, xx, yy, zz);
pass(1) = norm(imag(fVals(:)) - gVals(:) ) < 10*tol;
pass(2) = norm(imag(h) - g) < tol;

end
function pass = test_times( ) 
% Test times in diskfun 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f = diskfun(@(x,y) sin(x.*y)); 
pass(1) = norm(f.*f - f.^2, inf) < tol;

% Check for the case where f non-zero at the poles but g is zero at the
% poles.  The flag should indicate that the poles are zero.
z = diskfun(@(x,y) cos(x));
f = diskfun(@(x,y) sin(y));
g = z.*f;
pass(2) = g.nonZeroPoles == 0;
g = f.*z;
pass(3) = g.nonZeroPoles == 0;


end
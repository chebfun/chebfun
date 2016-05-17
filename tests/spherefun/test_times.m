function pass = test_times( ) 
% Test times in SPHEREFUN 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f = spherefun(@(x,y,z) sin(x.*y.*z)); 
pass(1) = norm(f.*f - f.^2, inf) < tol;

% Check for the case where f non-zero at the poles but g is zero at the
% poles.  The flag should indicate that the poles are zero.
z = spherefun(@(x,y,z) z);
f = spherefun(@(x,y,z) (1-z.^2));
g = z.*f;
pass(2) = g.nonZeroPoles == 0;
g = f.*z;
pass(3) = g.nonZeroPoles == 0;


end
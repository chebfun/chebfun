function pass = test_abs( ) 
% Test abs in SPHEREFUN 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f = spherefun(@(x,y,z) -(x.^2 + y.^2 + z.^2));
pass(1) = norm(abs(f) + f, inf) < tol;

end 
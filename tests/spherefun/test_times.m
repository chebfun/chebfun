function pass = test_times( ) 
% Test times in SPHEREFUN 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f = spherefun(@(x,y,z) sin(x.*y.*z)); 
pass(1) = norm(f.*f - f.^2, inf) < tol;

end
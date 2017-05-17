function pass = test_abs( ) 
% Test abs in DISKFUN 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f = diskfun(@(x,y) -(x.^2 + y.^2));
pass(1) = norm(abs(f) + f, inf) < tol;

end 
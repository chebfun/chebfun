function pass = test_median( ) 
% Test MEDIAN()

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

g = diskfun();
pass(1) = (isempty(median(g))); 

g = diskfun(@(x,y) 0*x+1); 
pass(2) = ( norm(median(g)-1) < tol);

end
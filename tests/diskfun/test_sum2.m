function pass = test_sum2( ) 
% Test diskfun sum2() command. 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f = @(x,y) 1 + x + y; 
g = diskfun(f);
exact_int = pi;
pass(1) = abs(sum2(g) - exact_int) < tol;

f = @(x,y) 10*exp(-x.*y) +12*x.^5.*cos(2*y)-10./(1+exp(x-.3));
g = diskfun(f);
exact_int = 14.16225689953454;
pass(2) = abs(sum2(g) - exact_int) < tol;

end 
function pass = test_integral2( ) 
% Test INTEGRAL2()

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

%test empty
pass(1) = (integral2(diskfun()) == 0); 

f = @(x,y) 1 + x + y; 
g = diskfun(f);
exact_int = pi;
pass(2) = abs(integral2(g) - exact_int) < tol;

f = @(x,y) 10*exp(-x.*y) +12*x.^5.*cos(2*y)-10./(1+exp(x-.3));
g = diskfun(f);
exact_int = 14.16225689953454;
pass(3) = abs(integral2(g) - exact_int) < tol;

end 


function pass = test_sum2( ) 
% Test spherefun sum2() command. 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f = @(x,y,z) 1 + x + y + z; 
g = spherefun(f);
exact_int = 4*pi;
pass(1) = abs(sum2(g) - exact_int) < tol;

f = @(x,y,z) 0.75*exp(-(9*x-2).^2/4 - (9*y-2).^2/4 - (9*z-2).^2/4)...
   + 0.75*exp(-(9*x+1).^2/49 - (9*y+1)/10 - (9*z+1)/10)...
   + 0.5*exp(-(9*x-7).^2/4 - (9*y-3).^2/4 - (9*z-5).^2/4)...
   - 0.2*exp(-(9*x-4).^2 - (9*y-7).^2 - (9*z-5).^2);
exact_int = 6.6961822200736179523;
g = spherefun(f);
pass(2) = abs(sum2(g) - exact_int) < tol;

end 
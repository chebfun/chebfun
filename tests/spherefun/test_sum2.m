function pass = test_sum2( ) 
% Test spherefun sum2() command. 

tol = 1e3*chebfunpref().techPrefs.eps;

f = @(x,y,z) 1 + x + y + z; 
g = spherefun( f ); 
exact_int = 4*pi;
pass(1) = abs( sum2( g ) - exact_int ) < tol;

end 
function pass = test_norm( ) 
% Test diskfun norm() command 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f = diskfun( @(x,y) 1+0*x );
pass(1) = abs( sqrt(sum2( f )) - norm( f ) ) < tol; 

f = diskfun( @(x,y) cos(x.*y) ); 
s = svd( f ); 
pass(2) = abs( sum(s.^2) - norm(f).^2 ) < tol; 

f = diskfun( @(x,y) x + y );  
pass(3) = abs( norm(f,inf) - sqrt(2) ) < tol; 

end
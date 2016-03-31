function pass = test_norm( ) 
% Test spherefun norm() command 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f = spherefun( @(x,y,z) 1+0*x );
pass(1) = abs( sqrt(sum2( f )) - norm( f ) ) < tol; 

f = spherefun( @(x,y,z) cos(x.*y.*z) ); 
s = svd( f ); 
pass(2) = abs( sum(s.^2) - norm(f).^2 ) < tol; 

f = spherefun( @(x,y,z) x + y + z );  
pass(3) = abs( norm(f,inf) - sqrt(3) ) < tol; 

end
function pass = test_helmholtzdecomp( ) 
% Test spherefunv/helmholtzdecomp code: 

tol = 1e5*chebfunpref().cheb2Prefs.chebfun2eps;

f = spherefunv( @(x,y,z) cos(x.*y.*z), @(x,y,z) sin(x+.1*y+5*z.^2),...
                                       @(x,y,z) x.*y.*z );                                 
f = tangent( f ); 
[u, v] = helmholtzdecomp( f ); 
pass(1) = norm( f - grad( u ) - curl( v )  ) < tol; 

f = spherefunv( @(x,y,z) sin((x-.1).*y.*z), @(x,y,z) cos(x+.50*y-1*z.^2),...
                                       @(x,y,z) -x.*y.^2.*z );  
f = tangent( f ); 
[u, v] = helmholtzdecomp( f ); 
pass(2) = norm( f - grad( u ) - curl( v )  ) < tol;

f = spherefunv( @(x,y,z) sin((x-.1).*y.*z), @(x,y,z) 1+0*x,...
                                       @(x,y,z) 0*x );  
f = tangent( f ); 
[u, v] = helmholtzdecomp( f ); 
pass(3) = norm( f - grad( u ) - curl( v )  ) < tol;

% Check empty spherefunv: 
[u, v] = helmholtzdecomp( spherefunv ); 
pass(4) = isempty( u ) && isempty( v );

end
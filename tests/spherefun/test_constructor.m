function pass = test_constructor( ) 
% Test the spherefun constructor 

% Get tolerance: 
tol = 2e3*chebfunpref().techPrefs.eps;

f = @(x,y,z) x.^2 + y.^2 + z.^2;
f = redefine_function_handle( f );
g = spherefun( f );
pass(1) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) exp(-cos(pi*(x+y+z)));
f = redefine_function_handle( f );
g = spherefun( f );
pass(2) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) 1-exp(x);
g = spherefun( f );
f = redefine_function_handle( f );
pass(3) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) exp(y);
g = spherefun( f );
f = redefine_function_handle( f );
pass(4) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) exp(z);
g = spherefun( f );
f = redefine_function_handle( f );
pass(5) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) cos(x.*y);
g = spherefun( f );
f = redefine_function_handle( f );
pass(6) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) sin(x.*y.*z);
g = spherefun( f );
f = redefine_function_handle( f );
pass(7) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) sin(x+ y.*z);
f = redefine_function_handle( f );
g = spherefun( f );
pass(8) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) sin(x+ y.*z) + 1;
f = redefine_function_handle( f );
g = spherefun( f );
pass(9) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) 0*x;
g = spherefun( f );
pass(10) = ( norm(g,inf) == 0 ); 

end
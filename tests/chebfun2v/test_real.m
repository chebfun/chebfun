function pass = test_real( pref ) 
% Test REAL

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 100*pref.eps; 

f = chebfun2v(@(x,y) cos(x.*y),@(x,y) cos(x.*y)); 
g = real( f ); 
pass(1) = ( norm( f - g ) < tol ); 

f = chebfun2v(@(x,y) cos(x.*y),@(x,y) cos(x.*y)); 
g = real( 1i*f ); 
pass(2) = ( norm( g ) < tol ); 

f1 = chebfun2v(@(x,y) cos(x.*y),@(x,y) cos(x.*y)); 
f2 = chebfun2v(@(x,y) sin(x + y.^2),@(x,y) sin(x + y.^2)); 
g = real( f1 + 1i*f2 ); 
pass(3) = ( norm( f1 - g ) < tol ); 

end
function pass = test_imag( pref ) 
% Test IMAG

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 10*pref.eps; 

f = chebfun2v(@(x,y) cos(x.*y),@(x,y) cos(x.*y)); 
g = imag( f ); 
pass(1) = ( norm( g ) < tol ); 

f = chebfun2v(@(x,y) cos(x.*y),@(x,y) cos(x.*y)); 
g = imag( 1i*f ); 
pass(2) = ( norm( f - g ) < 100*tol ); 

f1 = chebfun2v(@(x,y) cos(x.*y),@(x,y) cos(x.*y)); 
f2 = chebfun2v(@(x,y) sin(x + y.^2),@(x,y) sin(x + y.^2)); 
g = imag( f1 + 1i*f2 ); 
pass(3) = ( norm( f2 - g ) < 100*tol ); 

end
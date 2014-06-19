function pass = test_conj( pref ) 
% Test CONJ

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 100*pref.eps; 

f = chebfun2(@(x,y) cos(x.*y)); 
g = conj( f ); 
pass(1) = ( norm( f - g ) < tol ); 

f = chebfun2(@(x,y) cos(x.*y)); 
g = conj( 1i*f ); 
pass(2) = ( norm( 1i*f + g ) < 100*tol ); 

f1 = chebfun2(@(x,y) cos(x.*y)); 
f2 = chebfun2(@(x,y) sin(x + y.^2)); 
g = conj( f1 + 1i*f2 ); 
pass(3) = ( norm( f1 - 1i*f2 - g ) < tol ); 

end
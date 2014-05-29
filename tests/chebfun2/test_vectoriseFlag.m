function pass = test_vectoriseFlag( prefs ) 
% Test the vectorise flag in the constructor. 

if ( nargin < 1 ) 
    prefs = chebfunprefs(); 
end 

tol = prefs.cheb2Prefs.eps; 

% All these calls to the constructor should be the same: 
f1 = chebfun2(@(x,y) x); 
f2 = chebfun2(@(x,y) x, 'vectorize'); 
f3 = chebfun2(@(x,y) x, [-1,1,-1,1], 'vectorize'); 
f4 = chebfun2(@(x,y) x, 'vectorize', [-1,1,-1,1]); 
pass(1) = ( norm( f1 - f2 ) < tol ); 
pass(2) = ( norm( f1 - f3 ) < tol ); 
pass(3) = ( norm( f1 - f4 ) < tol ); 

% These should be the same: 
f1 = chebfun2(@(x,y) x.*y); 
f2 = chebfun2(@(x,y) x*y, 'vectorize'); 
pass(4) = ( norm( f1 - f2 ) < tol ); 

g = chebfun2(@(z) sum(z.^(0:9)),[-1 1 -1 1]*pi/2,'vectorise');
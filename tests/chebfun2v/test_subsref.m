function pass = test_subsref( pref )
% Test Chebfun2v subsref() command. 

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 1000*pref.eps; 


% Test recursive subsref:  
f = chebfun2(@(x, y) sin(x.*y));
F = [f ; f ; f];

G = G(1); 
exact = G.pivotValues; 

pass(1) = norm( exact - F(1).pivotValues ) < tol; 



end
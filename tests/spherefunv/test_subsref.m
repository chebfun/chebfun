function pass = test_subsref( pref )
% Test Spherefunv subsref() command. 

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 1000*pref.cheb2Prefs.chebfun2eps;


% Test recursive subsref:  
f = spherefun(@(x, y, z) sin(10*x.*y.*z));
F = [f ; f ; f];

G = F(1); 
exact = G.pivotValues; 

pass(1) = norm( exact - F(1).pivotValues ) < tol; 



end
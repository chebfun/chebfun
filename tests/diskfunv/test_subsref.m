function pass = test_subsref( pref )
% Test diskfunv subsref() command. 

if ( nargin == 0) 
    pref = chebfunpref; 
end

tol = 1000*pref.cheb2Prefs.chebfun2eps;

% Test recursive subsref:  
f = diskfun(@(x, y) sin(10*x.*y));
F = [f ; f ];

G = F(1); 
exact = G.pivotValues;  
pass(1) = norm( exact - F(1).pivotValues ) < tol; 



end
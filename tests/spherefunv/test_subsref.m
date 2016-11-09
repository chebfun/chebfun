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


F = spherefunv(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
g = chebfun3(@(x,y,z) x + y + z);
h_true = spherefun(@(x,y,z) x + y + z);
h = g(F);
pass(2) = ( norm(h - h_true) < tol );

G = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
H_true = spherefunv(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
H = G(F);
pass(3) = ( norm(H - H_true) < tol );

end
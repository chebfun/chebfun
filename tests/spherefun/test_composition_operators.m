function pass = test_composition_operators( pref )
% Check that composition operations are working. 

if ( nargin == 0 )
    pref = chebfunpref; 
end

tol = 1e3*pref.cheb2Prefs.chebfun2eps;
j = 1; 

f = spherefun(@(x,y,z) cos(x.*y) + sin(x.*y) + z -.1); 

% Multiplication
x = spherefun(@(x,y,z) x); 
y = spherefun(@(x,y,z) y);
z = spherefun(@(x,y,z) z);
exact = @(x,y,z) (cos(x.*y) + sin(x.*y) + z -.1).*sin(z.*(x-.1).*(y+.4)); 
g = spherefun(@(x,y,z) exact(x,y,z)); 
pass(j) = ( norm( g - f.*sin(z.*(x-.1).*(y+.4)) ) < tol ); j = j + 1;

% Cosine
exact = @(x,y,z) cos(cos(x.*y) + sin(x.*y) + z -.1); 
g = spherefun(@(x,y,z) exact(x,y,z)); 
pass(j) = ( norm( g - cos(f) ) < tol ); j = j + 1;

% Cosh
exact = @(x,y,z) cosh( cos(x.*y) + sin(x.*y) + z -.1); 
g = spherefun(@(x,y,z) exact(x,y,z)); 
pass(j) = ( norm( g - cosh(f) ) < tol ); j = j + 1;

% Sine
exact = @(x,y,z)  sin(cos(x.*y) + sin(x.*y) + z -.1); 
g = spherefun(@(x,y,z) exact(x,y,z)); 
pass(j) = ( norm( g - sin(f) ) < tol ); j = j + 1;

% Sinh
exact = @(x,y,z) sinh( cos(x.*y) + sin(x.*y) + z -.1); 
g = spherefun(@(x,y,z) exact(x,y,z)); 
pass(j) = ( norm( g - sinh(f) ) < tol ); j = j + 1;

% Multiple operations: 
f = spherefun(@(x,y,z) z + sin(pi*x.*y));
pass(j) = (norm(f+f+f-3*f) < 100*tol); j=j+1; 
pass(j) = (norm(f.*f-f.^2) < tol); j=j+1; 

end

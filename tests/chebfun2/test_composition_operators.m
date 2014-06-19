function pass = test_composition_operators( pref )
% Check that composition operations are working. 

if ( nargin == 0 )
    pref = chebfunpref; 
end

tol = 1000*pref.eps;
j = 1; 

f = chebfun2(@(x,y) cos(x.*y) + sin(x.*y) + y -.1); 

% Multiplication
x = chebfun2(@(x,y) x); 
y = chebfun2(@(x,y) y);
exact = @(x,y) (cos(x.*y) + sin(x.*y) + y -.1).*sin((x-.1).*(y+.4)); 
g = chebfun2(@(x,y) exact(x,y)); 
pass(j) = ( norm( g - f.*sin((x-.1).*(y+.4)) ) < tol ); j = j + 1;

% Cosine
exact = @(x,y) cos(cos(x.*y) + sin(x.*y) + y -.1); 
g = chebfun2(@(x,y) exact(x,y)); 
pass(j) = ( norm( g - cos(f) ) < tol ); j = j + 1;

% Cosh
exact = @(x,y) cosh( cos(x.*y) + sin(x.*y) + y -.1); 
g = chebfun2(@(x,y) exact(x,y)); 
pass(j) = ( norm( g - cosh(f) ) < tol ); j = j + 1;

% Sine
exact = @(x,y)  sin(cos(x.*y) + sin(x.*y) + y -.1); 
g = chebfun2(@(x,y) exact(x,y)); 
pass(j) = ( norm( g - sin(f) ) < tol ); j = j + 1;

% Sinh
exact = @(x,y) sinh( cos(x.*y) + sin(x.*y) + y -.1); 
g = chebfun2(@(x,y) exact(x,y)); 
pass(j) = ( norm( g - sinh(f) ) < tol ); j = j + 1;

% Multiple operations: 
f = chebfun2(@(x,y) sin(10*x.*y),[-1 2 -1 1]);
pass(j) = (norm(f+f+f-3*f) < 100*tol); j=j+1; 
pass(j) = (norm(f.*f-f.^2) < tol); j=j+1; 

end

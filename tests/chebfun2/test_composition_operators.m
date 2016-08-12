function pass = test_composition_operators( pref )
% Check that composition operations are working. 

if ( nargin == 0 )
    pref = chebfunpref; 
end

tol = 1000*pref.cheb2Prefs.chebfun2eps;
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

% g(f) with a CHEBFUN2 f and g a CHEBFUN
f = chebfun2(@(x,y) x+y);
g = chebfun(@(t) t.^2);
G = chebfun(@(t) [t, t.^2]);
h = compose(f, g);
h_expected = chebfun2(@(x,y) (x+y).^2);
pass(j) = ( norm(h - h_expected) < tol ); j=j+1;

f = chebfun2(@(x,y) x+y);
G = chebfun(@(t) [t, t.^2]);
H = compose(f, G);
H_expected = chebfun2v(@(x,y) x+y, @(x,y) (x+y).^2);
pass(j) = ( norm(H - H_expected) < tol ); j=j+1;

% g(f) with a CHEBFUN2 and g a CHEBMATRIX
f = chebfun2(@(x,y) x+y);
G = [chebfun(@(t) t); chebfun(@(t) t.^2)];
H = compose(f, G);
H_expected = chebfun2v(@(x,y) x+y, @(x,y) (x+y).^2);
pass(j) = ( norm(H - H_expected) < tol ); j=j+1;

end

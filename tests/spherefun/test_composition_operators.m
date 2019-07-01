function pass = test_composition_operators( pref )
% Check that composition operations are working. 

if ( nargin == 0 )
    pref = chebfunpref; 
end

tol = 1e3*pref.cheb2Prefs.chebfun2eps;

f = spherefun(@(x,y,z) cos(x.*y) + sin(x.*y) + z -.1); 

% Multiplication
x = spherefun(@(x,y,z) x); 
y = spherefun(@(x,y,z) y);
z = spherefun(@(x,y,z) z);
exact = @(x,y,z) (cos(x.*y) + sin(x.*y) + z -.1).*sin(z.*(x-.1).*(y+.4)); 
g = spherefun(@(x,y,z) exact(x,y,z)); 
pass(1) = ( norm( g - f.*sin(z.*(x-.1).*(y+.4)) ) < 100*tol );

% Cosine
exact = @(x,y,z) cos(cos(x.*y) + sin(x.*y) + z -.1); 
g = spherefun(@(x,y,z) exact(x,y,z)); 
pass(2) = ( norm( g - cos(f) ) < tol );

% Cosh
exact = @(x,y,z) cosh( cos(x.*y) + sin(x.*y) + z -.1); 
g = spherefun(@(x,y,z) exact(x,y,z)); 
pass(3) = ( norm( g - cosh(f) ) < tol );

% Sine
exact = @(x,y,z)  sin(cos(x.*y) + sin(x.*y) + z -.1); 
g = spherefun(@(x,y,z) exact(x,y,z)); 
pass(4) = ( norm( g - sin(f) ) < tol );

% Sinh
exact = @(x,y,z) sinh( cos(x.*y) + sin(x.*y) + z -.1); 
g = spherefun(@(x,y,z) exact(x,y,z)); 
pass(5) = ( norm( g - sinh(f) ) < tol );

% Multiple operations: 
f = spherefun(@(x,y,z) z + sin(pi*x.*y));
pass(6) = (norm(f+f+f-3*f) < 100*tol);
pass(7) = (norm(f.*f-f.^2) < tol);

% Composition of a spherefun with a chebfun (1 column):
f = spherefun(@(x,y,z) z + sin(pi*x.*y));
g = chebfun(@(t) t.^2, [ -1.5, 1.5 ]);
h_true = spherefun(@(x,y,z) (z + sin(pi*x.*y)).^2);
h = compose(f, g);
pass(8) = ( norm(h - h_true) < tol );

% Composition of a spherefun with a chebfun (3 columns):
f = spherefun(@(x,y,z) z + sin(pi*x.*y));
G = chebfun(@(t) [ t.^2, t, -t.^2 ], [ -1.5, 1.5 ]);
H_true = spherefunv(@(x,y,z) (z + sin(pi*x.*y)).^2, ...
    @(x,y,z) z + sin(pi*x.*y), @(x,y,z) -(z + sin(pi*x.*y)).^2);
H = compose(f, G);
pass(9) = ( norm(H - H_true) < tol );

% Composition of a spherefun with a Chebfun2:
f = spherefun(@(x,y,z) z + sin(pi*x.*y));
g = chebfun2(@(x,y) x + y, [-2, 2, -2, 2]);
h_true = f;
h = compose(f, g);
pass(10) = ( norm(h - h_true) < tol );

% ... and with a Chebfun2v:
G = [g; g; g];
H_true = [f; f; f];
H = compose(f, G);
pass(11) = ( norm(H - H_true) < tol );

end

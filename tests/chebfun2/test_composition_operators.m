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
pass(1) = ( norm( g - f.*sin((x-.1).*(y+.4)) ) < tol );

% Cosine
exact = @(x,y) cos(cos(x.*y) + sin(x.*y) + y -.1); 
g = chebfun2(@(x,y) exact(x,y)); 
pass(2) = ( norm( g - cos(f) ) < tol );

% Cosh
exact = @(x,y) cosh( cos(x.*y) + sin(x.*y) + y -.1); 
g = chebfun2(@(x,y) exact(x,y)); 
pass(3) = ( norm( g - cosh(f) ) < tol );

% Sine
exact = @(x,y)  sin(cos(x.*y) + sin(x.*y) + y -.1); 
g = chebfun2(@(x,y) exact(x,y)); 
pass(4) = ( norm( g - sin(f) ) < tol );

% Sinh
exact = @(x,y) sinh( cos(x.*y) + sin(x.*y) + y -.1); 
g = chebfun2(@(x,y) exact(x,y)); 
pass(5) = ( norm( g - sinh(f) ) < tol );

% Multiple operations: 
f = chebfun2(@(x,y) sin(10*x.*y),[-1 2 -1 1]);
pass(6) = (norm(f+f+f-3*f) < 100*tol);
pass(7) = (norm(f.*f-f.^2) < tol);

% Compose a CHEBFUN2 f with a CHEBFUN g (one column):
f = chebfun2(@(x,y) x + y);
g = chebfun(@(t) t.^2, [-2, 2]);
h = compose(f, g);
h_true = chebfun2(@(x,y) (x + y).^2);
pass(8) = ( norm(h - h_true) < tol );
pass(9) = ~isPeriodicTech(h);

% Compose a CHEBFUN2 f with a CHEBFUN G (two columns):
f = chebfun2(@(x,y) x + y);
G = chebfun(@(t) [t, t.^2], [-2, 2]);
H = compose(f, G);
H_true = chebfun2v(@(x,y) x + y, @(x,y) (x + y).^2);
pass(10) = ( norm(H - H_true) < tol );

% Compose a periodic CHEBFUN2 with a CHEBFUN g:
f = chebfun2(@(x,y) cos(pi*x) .* sin(y), [ -1, 1, -pi, pi ], 'trig');
g = chebfun(@(t) t.^2);
h = compose(f, g);
pass(11) = isPeriodicTech(h);
pass(12) = ( norm(h.domain - f.domain) < tol );

g = chebfun(@(t) [ t, cos(t) ]);
h = compose(f, g);
pass(13) = isPeriodicTech(h);

% Compose a complex-valued CHEBFUN2 f with a CHEBFUN2 g:
f = chebfun2(@(x,y) x + 1i*y);
g = chebfun2(@(x,y) x.^2 + y.^2);
h_true = chebfun2(@(x,y) x.^2 + y.^2);
h = compose(f, g);
pass(14) = ( norm(h - h_true) < tol );

% Compose a periodic complex CHEBFUN2 f with a CHEBFUN2 g:
f = chebfun2(@(x,y) exp(1i * pi * x), 'trig');
g = chebfun2(@(x,y) x);
h = compose(f, g);
pass(15) = isPeriodicTech(h);
G = chebfun2v(@(x,y) x, @(x,y) y);
H = compose(f, G);
pass(16) = isPeriodicTech(H);

end

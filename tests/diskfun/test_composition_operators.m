function pass = test_composition_operators( pref )
% Check that composition operations are working.

if ( nargin == 0 )
    pref = chebfunpref;
end

tol = 1e3*pref.cheb2Prefs.chebfun2eps;

f = diskfun(@(x,y) cos(x.*y) + sin(x.*(y)) );

% Multiplication
x = diskfun(@(x,y) x);
y = diskfun(@(x,y) y);
exact = @(x,y) (cos(x.*y) + sin(x.*y)).*sin((x-.1).*(y+.4));
g = diskfun(@(x,y) exact(x,y));
pass(1) = ( norm( g - f.*sin((x-.1).*(y+.4)) ) < tol );


% Cosine
exact = @(x,y) cos(cos(x.*y) + sin(x.*y));
g = diskfun(@(x,y) exact(x,y));
pass(2) = ( norm( g - cos(f) ) < tol );

% Cosh
exact = @(x,y) cosh( cos(x.*y) + sin(x.*y) );
g = diskfun(@(x,y) exact(x,y));
pass(3) = ( norm( g - cosh(f) ) < tol );

% Sine
exact = @(x,y)  sin(cos(x.*y) + sin(x.*y));
g = diskfun(@(x,y,z) exact(x,y));
pass(4) = ( norm( g - sin(f) ) < tol );

% Sinh
exact = @(x,y) sinh( cos(x.*y) + sin(x.*y));
g = diskfun(@(x,y) exact(x,y));
pass(5) = ( norm( g - sinh(f) ) < tol );

% Multiple operations:
f = diskfun(@(x,y)  sin(pi*x.*y));
pass(6) = (norm(f+f+f-3*f) < 100*tol);
pass(7) = (norm(f.*f-f.^2) < tol);

% Composition with a CHEBFUN with one column:
f = diskfun(@(x,y) x + y);
g = chebfun(@(t) t.^2, [-2, 2]);
h = compose(f, g);
h_true = diskfun(@(x,y) (x + y).^2);
pass(8) = ( norm(h - h_true) < tol );

% ... and a CHEBFUN with two columns:
G = chebfun(@(t) [t.^2, exp(t)], [-2, 2]);
H = compose(f, G);
H_true = diskfunv(@(x,y) (x + y).^2, @(x,y) exp(x + y));
pass(9) = ( norm(H - H_true) < tol );

% Composition with a CHEBFUN2 and a CHEBFUN2V:
f = diskfun(@(x,y) x + y);
g = chebfun2(@(x,y) x.^2 + y.^2, [-2, 2, -2, 2]);
h = compose(f, g);
h_true = diskfun(@(x,y) (x + y).^2);
pass(10) = ( norm(h - h_true) < tol );

H = compose(f, [g; g]);
H_true = [h_true; h_true];
pass(11) = ( norm(H - H_true) < tol );

end

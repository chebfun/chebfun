function pass = test_compose(pref)
% Check that composition operations are working. 

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1000*pref.cheb2Prefs.chebfun2eps;

%% F = CHEBFUN2V, g = CHEBFUN2
F = chebfun2v(@(x,y) x, @(x,y) y, [ 0, 1, 0, 1 ]);
g = chebfun2(@(x,y) x + y);
h = compose(F, g);
h_true = chebfun2(@(x,y) x + y, [ 0, 1, 0, 1 ]);
pass(1) = ( norm(h - h_true) < tol );

%% F = CHEBFUN2V, G = CHEBFUN2V with 2 components
F = chebfun2v(@(x,y) x, @(x,y) y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y);
H = compose(F, G);
pass(2) = ( norm(H - G) < tol );

%% F = CHEBFUN2V, G = CHEBFUN2V with 3 components
F = chebfun2v(@(x,y) x, @(x,y) y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y, @(x,y) x .* y);
H = compose(F, G);
pass(3) = ( norm(H - G) < tol );

%% F = CHEBFUN2V, g = CHEBFUN3
F = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y);
g = chebfun3(@(x,y,z) x + y + z, [ -1, 1, -1, 1, -2, 2 ]);
h = compose(F, g);
h_true = chebfun2(@(x,y) 2*x + 2*y);
pass(4) = ( norm(h - h_true) < tol );

%% F = CHEBFUN2V, G = CHEBFUN3V with 2 components
F = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y);
G = chebfun3v(@(x,y,z) x + y, @(x,y,z) z - x, [ -1, 1, -1, 1, -2, 2 ]);
H = compose(F, G);
H_true = chebfun2v(@(x,y) x+y, @(x,y) y);
pass(5) = ( norm(H - H_true) < tol );

%% F = CHEBFUN2V, G = CHEBFUN3V with 3 components
F = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y);
G = chebfun3v(@(x,y,z) x + y, @(x,y,z) x - y, @(x,y,z) z, [ -1, 1, -1, 1, -2, 2 ]);
H = compose(F, G);
H_true = chebfun2v(@(x,y) x+y, @(x,y) x - y, @(x,y) x + y);
pass(6) = ( norm(H - H_true) < tol );

%% F = peridioc CHEBFUN2V, G = CHEBFUN2 and CHEBFUN2V
f1 = chebfun2(@(x,y) cos(pi*x) .* sin(pi*y), 'trig');
f2 = chebfun2(@(x,y) cos(pi*x) .* cos(pi*y), 'trig');
F = [f1; f2];
g = chebfun2(@(x,y) x.^2 + y.^2);
h = g(F);
pass(7) = isPeriodicTech(h);

G = chebfun2v(@(x,y) x.^2 + y.^2, @(x,y) x);
H = G(F);
pass(8) = isPeriodicTech(H);

%% F = periodic CHEBFUN2V and G = CHEBFUN3 or CHEBFUN3V
F = [f1; f2; f1];
g = chebfun3(@(x,y,z) exp(x) .* cos(y) + 1);
h = g(F);
pass(9) = isPeriodicTech(h);
G = chebfun3v(@(x,y,z) exp(x), @(x,y,z) x.^2 + y.^2, @(x,y,z) z-2);
H = G(F);
pass(10) = isPeriodicTech(H);

end
function pass = test_compose(pref)
% Check that composition operations are working. 

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1000*pref.cheb2Prefs.chebfun2eps;
j = 1;

%% F = CHEBFUN2V, g = CHEBFUN2
F = chebfun2v(@(x,y) x, @(x,y) y, [ 0, 1, 0, 1 ]);
g = chebfun2(@(x,y) x + y);
h = compose(F, g);
h_true = chebfun2(@(x,y) x + y, [ 0, 1, 0, 1 ]);
pass(j) = ( norm(h - h_true) < tol );
j = j+1;

%% F = CHEBFUN2V, G = CHEBFUN2V with 2 components
F = chebfun2v(@(x,y) x, @(x,y) y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y);
H = compose(F, G);
pass(j) = ( norm(H - G) < tol );
j=j+1;

%% F = CHEBFUN2V, G = CHEBFUN2V with 3 components
F = chebfun2v(@(x,y) x, @(x,y) y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y, @(x,y) x .* y);
H = compose(F, G);
pass(j) = ( norm(H - G) < tol );
j=j+1;

%% F = CHEBFUN2V, g = CHEBFUN3
F = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y);
g = chebfun3(@(x,y,z) x + y + z);
h = compose(F, g);
h_true = chebfun2(@(x,y) 2*x + 2*y);
pass(j) = ( norm(h - h_true) < tol );
j = j+1;

%% F = CHEBFUN2V, G = CHEBFUN3V with 2 components
F = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y);
G = chebfun3v(@(x,y,z) x + y, @(x,y,z) z - x);
H = compose(F, G);
H_true = chebfun2v(@(x,y) x+y, @(x,y) y);
pass(j) = ( norm(H - H_true) < tol );
j = j+1;

%% F = CHEBFUN2V, G = CHEBFUN3V with 3 components
F = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y);
G = chebfun3v(@(x,y,z) x + y, @(x,y,z) x - y, @(x,y,z) z);
H = compose(F, G);
H_true = chebfun2v(@(x,y) x+y, @(x,y) x - y, @(x,y) x + y);
pass(j) = ( norm(H - H_true) < tol );
j = j+1;

%% F = peridioc CHEBFUN2V, G = CHEBFUN2 and CHEBFUN2V
f1 = chebfun2(@(x,y) cos(pi*x) .* sin(pi*y), 'trig');
f2 = chebfun2(@(x,y) cos(pi*x) .* cos(pi*y), 'trig');
F = [f1; f2];
g = chebfun2(@(x,y) x.^2 + y.^2);
h = g(F);
pass(j) = isPeriodicTech(h);
j = j + 1;

G = chebfun2v(@(x,y) x.^2 + y.^2, @(x,y) x);
H = G(F);
pass(j) = isPeriodicTech(H);
j = j + 1;

%% F = periodic CHEBFUN2V and G = CHEBFUN3 or CHEBFUN3V
F = [f1; f2; f1];
g = chebfun3(@(x,y,z) exp(x) .* cos(y) + 1);
h = g(F);
pass(j) = isPeriodicTech(h);
j = j + 1;
G = chebfun3v(@(x,y,z) exp(x), @(x,y,z) x.^2 + y.^2, @(x,y,z) z-2);
H = G(F);
pass(j) = isPeriodicTech(H);

end
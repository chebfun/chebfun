function pass = test_compose(pref)
% Check that composition operations are working.

if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;
j = 1;

%% Compose a CHEBFUN3V F (3 components) with a CHEBFUN3 g:
F = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z, [ 0, 1, 0, 1, 0, 1 ]);
g = chebfun3(@(x,y,z) x + y + z);
h = compose(F, g);
h_true = chebfun3(@(x,y,z) x + y + z, [ 0, 1, 0, 1, 0, 1 ]);
pass(j) = ( norm(h - h_true) < tol );
j = j+1;

%% Compose a CHEBFUN3V F (3 components) with a CHEBFUN3V G (2 components):
F = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
G = chebfun3v(@(x,y,z) x + y + z, @(x,y,z) x + z);
H = compose(F, G);
pass(j) = ( norm(H - G) < tol );
j = j+1;

%% Compose a CHEBFUN3V F (3 components) with a CHEBFUN3V G (3 components):
F = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
G = chebfun3v(@(x,y,z) x + y + z, @(x,y,z) x + z, @(x,y,z) x - y);
H = compose(F, G);
pass(j) = ( norm(H - G) < tol );
j = j+1;

%% Compose a CHEBFUN3V F (2 components) with a CHEBFUN2 g:
F = chebfun3v(@(x,y,z) x, @(x,y,z) y);
g = chebfun2(@(x,y) x + y);
h = compose(F, g);
h_true = chebfun3(@(x,y,z) x + y);
pass(j) = ( norm(h - h_true) < tol );
j = j+1;

%% Compose a CHEBFUN3V F (2 components) with a CHEBFUN2V G (2 components):
F = chebfun3v(@(x,y,z) x, @(x,y,z) y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y);
H = compose(F, G);
H_true = chebfun3v(@(x,y,z) x + y, @(x,y,z) x - y);
pass(j) = ( norm(H - H_true) < tol );
j = j+1;

%% Compose a CHEBFUN3V F (2 components) with a CHEBFUN2V G (3 components):
F = chebfun3v(@(x,y,z) x, @(x,y,z) y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y, @(x,y) x .* y);
H = compose(F, G);
H_true = chebfun3v(@(x,y,z) x + y, @(x,y,z) x - y, @(x,y,z) x .* y);
pass(j) = ( norm(H - H_true) < tol );
j = j + 1;

%% Compose a periodic CHEBFUN3V F (2 components) with a CHEBFUN2 and CHEBFUN2V:
f1 = chebfun3(@(x,y,z) cos(pi*x).*sin(pi*y), 'trig');
f2 = chebfun3(@(x,y,z) sin(pi*x).*sin(pi*y), 'trig');
F = [ f1; f2 ];
g = chebfun2(@(x,y) x, [ -1, 1, -1, 1 ]);
h = compose(F, g);
pass(j) = isPeriodicTech(h);
j = j + 1;
G = chebfun2v(@(x,y) x, @(x,y) y, [ -1, 1, -1, 1 ]);
H = compose(F, G);
pass(j) = isPeriodicTech(H);
j = j + 1;

%% Compose a periodic CHEBFUN3V F (3 components) with a CHEBFUN3 and CHEBFUN3V:
f1 = chebfun3(@(x,y,z) cos(pi*x).*sin(pi*y), 'trig');
f2 = chebfun3(@(x,y,z) sin(pi*x).*sin(pi*y), 'trig');
F = [ f1; f2; f1 ];
g = chebfun3(@(x,y,z) x, [ -1, 1, -1, 1, -1, 1 ]);
h = compose(F, g);
pass(j) = isPeriodicTech(h);
j = j + 1;
G = chebfun3v(@(x,y,z) x, @(x,y,z) y, [ -1, 1, -1, 1, -1, 1 ]);
H = compose(F, G);
pass(j) = isPeriodicTech(H);

end
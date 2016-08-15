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

end
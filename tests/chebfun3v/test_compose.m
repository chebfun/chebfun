function pass = test_compose(pref)
% Check that composition operations are working.

if ( nargin == 0 )
    pref = chebfunpref;
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;
j = 1;

%% F = CHEBFUN3V, g = CHEBFUN3
F = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
g = chebfun3(@(x,y,z) x + y + z);
h = compose(F, g);
pass(j) = ( norm(h - g) < tol );
j = j+1;

%% F = CHEBFUN3V, G = CHEBFUN3V with 2 components
F = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
G = chebfun3v(@(x,y,z) x + y + z, @(x,y,z) x + z);
H = compose(F, G);
pass(j) = ( norm(H - G) < tol );
j = j+1;

%% F = CHEBFUN3V, G = CHEBFUN3V with 3 components
F = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
G = chebfun3v(@(x,y,z) x + y + z, @(x,y,z) x + z, @(x,y,z) x - y);
H = compose(F, G);
pass(j) = ( norm(H - G) < tol );
j = j+1;

%% F = CHEBFUN3V, g = CHEBFUN2
F = chebfun3v(@(x,y,z) x, @(x,y,z) y);
g = chebfun2(@(x,y) x + y);
h = compose(F, g);
h_expected = chebfun3(@(x,y,z) x + y);
pass(j) = ( norm(h - h_expected) < tol );
j = j+1;

%% F = CHEBFUN3V, G = CHEBFUN2V with 2 components
F = chebfun3v(@(x,y,z) x, @(x,y,z) y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y);
H = compose(F, G);
H_expected = chebfun3v(@(x,y,z) x + y, @(x,y,z) x - y);
pass(j) = ( norm(H - H_expected) < tol );
j = j+1;

%% F = CHEBFUN3V, G = CHEBFUN2V with 3 components
F = chebfun3v(@(x,y,z) x, @(x,y,z) y);
G = chebfun2v(@(x,y) x + y, @(x,y) x - y, @(x,y) x .* y);
H = compose(F, G);
H_expected = chebfun3v(@(x,y,z) x + y, @(x,y,z) x - y, @(x,y,z) x .* y);
pass(j) = ( norm(H - H_expected) < tol );
j = j+1;

end
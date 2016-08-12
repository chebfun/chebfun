function pass = test_compose(pref)
% Check that composition operations are working. 

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;


%% F = CHEBFUN3V, g = CHEBFUN3
F = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
g = chebfun3(@(x,y,z) x + y + z);
h = compose(F, g);
pass(1) = ( norm(h - g) < tol );

%% F = CHEBFUN3V, G = CHEBFUN3V
F = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
G = chebfun3v(@(x,y,z) x + y + z, @(x,y,z) x+z);
H = compose(F, G);
pass(2) = ( norm(H - G) < tol );

%% F = CHEBFUN3V, g = CHEBFUN2
F = chebfun3v(@(x,y,z) x, @(x,y,z) y);
g = chebfun2(@(x,y) x + y);
h = compose(F, g);
expected = chebfun3(@(x,y,z) x+y);
pass(3) = ( norm(h - expected) < tol );

%% F = CHEBFUN3V, G = CHEBFUN2V
F = chebfun3v(@(x,y,z) x, @(x,y,z) y);
G = chebfun2v(@(x,y) x + y, @(x,y) x-y);
H = compose(F, G);
expected = chebfun3v(@(x,y,z) x+y, @(x,y,z) x-y);
pass(4) = ( norm(H - expected) < tol );

end
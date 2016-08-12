function pass = test_compose(pref)
% Check that composition operations are working. 

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1000*pref.cheb2Prefs.chebfun2eps;
j = 1;

%% F = CHEBFUN2V, g = CHEBFUN2
F = chebfun2v(@(x,y) x, @(x,y) y);
g = chebfun2(@(x,y) x + y);
h = compose(F, g);
pass(j) = ( norm(h - g) < tol );
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
h_expected = chebfun2(@(x,y) 2*x + 2*y);
pass(j) = ( norm(h - h_expected) < tol );
j = j+1;

%% F = CHEBFUN2V, G = CHEBFUN3V with 2 components
F = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y);
G = chebfun3v(@(x,y,z) x + y, @(x,y,z) z - x);
H = compose(F, G);
H_expected = chebfun2v(@(x,y) x+y, @(x,y) y);
pass(j) = ( norm(H - H_expected) < tol );
j = j+1;

%% F = CHEBFUN2V, G = CHEBFUN3V with 3 components
F = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y);
G = chebfun3v(@(x,y,z) x + y, @(x,y,z) x - y, @(x,y,z) z);
H = compose(F, G);
H_expected = chebfun2v(@(x,y) x+y, @(x,y) x - y, @(x,y) x + y);
pass(j) = ( norm(H - H_expected) < tol );
j = j+1;

end
function pass = test_compose(pref)
% Test COMPOSE().

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1e3*pref.cheb2Prefs.chebfun2eps;

F = diskfunv(@(x,y) x, @(x,y) y);
g = chebfun2(@(x,y) x + y);
h = compose(F, g);
h_true = diskfun(@(x,y) x + y);
pass(1) = ( norm(h - h_true) < tol );

G = chebfun2v(@(x,y) x + y, @(x,y) x - y);
H = compose(F, G);
H_true = diskfunv(@(x,y) x + y, @(x,y) x - y);
pass(2) = ( norm(H - H_true) < tol );

end

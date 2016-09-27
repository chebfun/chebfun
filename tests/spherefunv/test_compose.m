function pass = test_compose(pref)
% Test COMPOSE().

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1e3*pref.cheb2Prefs.chebfun2eps;

F = spherefunv(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
g = chebfun3(@(x,y,z) x + y + z);
h_true = spherefun(@(x,y,z) x + y + z);
h = compose(F, g);
pass(1) = ( norm(h - h_true) < tol );

G = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
H_true = spherefunv(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
H = compose(F, G);
pass(2) = ( norm(H - H_true) < tol );

end

function pass = test_compose(pref)
% Test COMPOSE().

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1e3*pref.cheb2Prefs.chebfun2eps;
j = 1; 

F = spherefunv(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
g = chebfun3(@(x,y,z) x + y + z);
h_true = spherefun(@(x,y,z) x + y + z);
h = compose(F, g);
pass(j) = ( norm(h - h_true) < tol );
j = j + 1;

G = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
H_true = spherefunv(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
H = compose(F, G);
pass(j) = ( norm(H - H_true) < tol );

end

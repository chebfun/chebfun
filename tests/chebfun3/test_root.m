function pass = test_root(pref)
% Check chebfun3/root.

if ( nargin < 1 ) 
    pref = chebfunpref;
end
tol = 10*pref.cheb3Prefs.chebfun3eps;

% Try a rootfinding problem from Chebfun3 guide:
f = chebfun3(@(x,y,z) y-x.^2);
g = chebfun3(@(x,y,z) z-x.^3);
h = chebfun3(@(x,y,z) cos(exp(x.*sin(-2+y+z))));
r = root(f, g, h);

%%
% Is the root accurate?
pass(1) = f(r(1), r(2), r(3)) < tol;
pass(2) = g(r(1), r(2), r(3)) < tol;
pass(3) = h(r(1), r(2), r(3)) < tol;


%% Another trivial test. See #1900:
tol = 1e-2;
f = chebfun3(@(x,y,z) x.^2-tol);
g = chebfun3(@(x,y,z) y.^2-tol);
h = chebfun3(@(x,y,z) z.^2-tol);
r = root(f, g, h);
pass(4) = f(r(1), r(2), r(3)) < tol;
pass(5) = g(r(1), r(2), r(3)) < tol;
pass(6) = h(r(1), r(2), r(3)) < tol;
end
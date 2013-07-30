function pass = test_merge(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

tol = 10*pref.chebfun.eps;

% Test something easy (the example from docs):
pref2 = pref;
pref2.chebfun.splitting = 1;
f = chebfun(@(x) abs(x), pref2);
f2 = f.^2;
g = merge(f2);
pass(1) = numel(g.funs) == 1;
xx = linspace(-1, 1, 100);
pass(2) = norm(f2(xx) - g(xx), 'inf') < tol;

% Test selective merge on many points:
pref2 = pref;
pref2.chebfun.splitting = 1;
f = chebfun(@(x) sin(10*pi*x), [-1:.5:2], pref2);
[g, mergedPts] = merge(f, [2 4:6]);
pass(3) = numel(g.funs) == 2;
xx = linspace(-1, 2, 100);
pass(4) = norm(f(xx) - g(xx), 'inf') < 10*tol;
pass(5) = all(mergedPts == [2 4:6]);

% Test a non-smooth function:
x = chebfun('x', [-1 0 1], pref);
f = sin(x);
g = abs(x-.5);
h = f + g;
[p, mergedPts] = merge(h);
pass(6) = numel(p.funs) == 2 && norm(p.domain - [-1 .5 1], inf) < tol;
pass(7) = norm(h(xx) - p(xx), 'inf') < 10*tol;
pass(8) = all(mergedPts == 2);

% Test a function with an impulse:
x = chebfun('x', [-1 0 1], pref);
x.impulses = [-1 0 1 ; 0 1 0];
pass(9) = isequal(merge(x), x);

end
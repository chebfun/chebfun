% Test file for @chebfun/merge.m.

function pass = test_merge(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

seedRNG(6178);

% Test something easy (the example from docs):
pref2 = pref;
pref2.splitting = 1;
f = chebfun(@(x) x.^2, [-1 0 1], pref2);
g = merge(f);
pass(1) = all(size(g.funs) == 1);
xx = linspace(-1, 1, 100);
pass(2) = norm(feval(f, xx) - feval(g, xx), 'inf') < 10*eps;

% Test selective merge on many points:
pref2 = pref;
pref2.splitting = 1;
f = chebfun(@(x) sin(10*pi*x), [-1:.5:2], pref2);
[g, mergedPts] = merge(f, [2 4:6]);
pass(3) = all(size(g.funs) == [1,2]);
xx = linspace(-1, 2, 100);
err = norm(feval(f, xx) - feval(g, xx), 'inf');
tol = 1e2*eps;

pass(4) = err < tol;
pass(5) = all(mergedPts == [2 4:6]);

% Test a non-smooth function:
h = chebfun(@(x) sin(x) + abs(x - 0.5), [-1 0 1], 'splitting', 'on');
[p, mergedPts] = merge(h);
pass(6) = numel(p.funs) == 2 && norm(p.domain - [-1 .5 1], inf) < 10*eps;
pass(7) = norm(feval(h, xx) - feval(p, xx), 'inf') < 10*eps;
pass(8) = all(mergedPts == 2);

% Test an array-valued CHEBFUN:
f = chebfun(@(x) [x, x.^2], [-1 -.1 -.1+eps 0 1]);
g = merge(f);
pass(9) = numel(g.domain) == 2 && all(g.domain == [-1 1]);

% Test row CHEBFUNs:
f = chebfun(@(x) x, [-1 0 1]);
g = merge(f.');
pass(10) = isequal(g.domain, [-1 1]) && g.isTransposed;

f = chebfun(@(x) abs(x), [-1 0 1]);
g = merge(f.');
pass(11) = isequal(g.domain, [-1 0 1]) && g.isTransposed;

%% Test for singular function:
% Set the domain:
dom = [-2 7];

pow1 = -1;
pow2 = -1;
op = @(x) (x - dom(1)).^pow1.*sin(10*x).*(x - dom(2)).^pow2;
pref = chebfunpref();
f = chebfun(op, dom, 'exps', [-1 -1]);
g = addBreaksAtRoots(f);
h = merge(g);

% Checking domain:
domCheck = [dom(1)+0.1 dom(2)-0.1];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

vals_h = feval(h, x);
vals_exact = feval(op, x);
err = vals_h - vals_exact;
pass(12) = (norm(err, inf) < 2e4*get(h, 'vscale')*eps);


%% Test for function defined on unbounded domain:

dom = [0:10:100 Inf];
domCheck = [0 100];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) 0.75+sin(10*x)./exp(x);
f = chebfun(op, dom, 'splitting', 'on');
g = merge(f);
gVals = feval(g, x);
gExact = op(x);
err = norm(gVals - gExact, inf);
tol = 1e3*eps*vscale(f);
pass(13) = err < tol;

%% Test for vertial scale invariance:
g = chebfun(@(x) 1./(1+25*x.^2) + 1e-14*sign(x), [-1 0 1]);
mg = merge(g);
h = chebfun(@(x) 100 * (1./(1+25*x.^2) + 1e-14*sign(x)), [-1 0 1]);
mh = merge(h);
pass(14) = numel(mg.funs) == numel(mh.funs);

end

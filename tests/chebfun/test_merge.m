% Test file for @chebfun/merge.m.

function pass = test_merge(pref)

if ( nargin == 0 )
    pref = chebpref();
end

% Test something easy (the example from docs):
pref2 = pref;
pref2.enableBreakpointDetection = 1;
f = chebfun(@(x) x.^2, [-1 0 1], pref2);
g = merge(f);
pass(1) = all(size(g.funs) == 1);
xx = linspace(-1, 1, 100);
pass(2) = norm(feval(f, xx) - feval(g, xx), 'inf') < epslevel(f);

% Test selective merge on many points:
pref2 = pref;
pref2.enableBreakpointDetection = 1;
f = chebfun(@(x) sin(10*pi*x), [-1:.5:2], pref2);
[g, mergedPts] = merge(f, [2 4:6]);
pass(3) = all(size(g.funs) == [1,2]);
xx = linspace(-1, 2, 100);
pass(4) = norm(feval(f, xx) - feval(g, xx), 'inf') < 10*epslevel(f);
pass(5) = all(mergedPts == [2 4:6]);

% Test a non-smooth function:
h = chebfun(@(x) sin(x) + abs(x - 0.5), [-1 0 1], 'splitting', 'on');
[p, mergedPts] = merge(h);
pass(6) = numel(p.funs) == 2 && norm(p.domain - [-1 .5 1], inf) < 10*eps;
pass(7) = norm(feval(h, xx) - feval(p, xx), 'inf') < 10*epslevel(h);
pass(8) = all(mergedPts == 2);

% Test a function with an impulse:
x = chebfun('x', [-1 0 1], pref);
% Artificially insert an impulse:
x.impulses = cat(3, [-1 0 1].', [0 1 0].');
% The break at zero should not be removed:
pass(9) = isequal(merge(x), x);

% Test an array-valued CHEBFUN:
f = chebfun(@(x) [x, x.^2], [-1 -.1 -.1+eps 0 1]);
g = merge(f);
pass(10) = numel(g.domain) == 2 && all(g.domain == [-1 1]);

%% Test for singular function:
% Set the domain:
dom = [-2 7];

pow1 = -1;
pow2 = -1;
op = @(x) (x - dom(1)).^pow1.*sin(10*x).*(x - dom(2)).^pow2;
pref = chebpref();
pref.singPrefs.exponents = [-1 -1];
f = chebfun(op, dom, pref);
g = addBreaksAtRoots(f);
h = merge(g);

% Checking domain:
domCheck = [dom(1)+0.1 dom(2)-0.1];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

vals_h = feval(h, x);
vals_exact = feval(op, x);
err = vals_h - vals_exact;
pass(11) = (norm(err, inf) < 5e1*get(h, 'vscale')*get(h, 'epslevel'));

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
err = gVals - gExact;
pass(12) = norm(err, inf) < epslevel(f)*vscale(f);

end
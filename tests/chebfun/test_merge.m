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

% Test an array-valued CHEBFUN:
f = chebfun(@(x) [x, x.^2], [-1 -.1 -.1+eps 0 1]);
g = merge(f);
pass(9) = numel(g.domain) == 2 && all(g.domain == [-1 1]);

end

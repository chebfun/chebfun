% Test file for @chebfun/not.m.

function pass = test_not(pref)

if (nargin < 1)
    pref = chebpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Check a few basic cases.
f = chebfun(@(x) sin(2*x), [-1 1], pref);
g = not(f);
pass(1) = isequal(g.impulses, [0 1 0].') && all(feval(g, x) == 0);

f = chebfun(@(x) exp(x), [-1 -0.5 0.5 1], pref);
g = not(f);
pass(2) = all(g.impulses == 0) && all(feval(g, x) == 0);

f = chebfun(@(x) 0*x, [-1 -0.5 0.5 1], pref);
g = not(f);
pass(3) = all(g.impulses == 1) && all(feval(g, x) == 1);

f = chebfun(@(x) 0*x + 1, [-1 -0.5 0.5 1], pref);
g = not(f);
pass(4) = all(g.impulses == 0) && all(feval(g, x) == 0);

% Check complex chebfun.
f = chebfun(@(x) exp(1i*x).*sin(x), [-1 -0.5 0.5 1], pref);
g = not(f);
pass(5) = isequal(g.impulses, [0 0 1 0 0].') && all(feval(g, x) == 0);

% Check array-valued chebfun.
f = chebfun(@(x) [sin(x) exp(x)], [-1 -0.5 0.5 1], pref);
g = not(f);
pass(6) = isequal(g.impulses, [0 0 1 0 0 ; 0 0 0 0 0].') && ...
    all(all(feval(g, x) == 0));

end

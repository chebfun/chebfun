% Test file for @chebfun/or.m.

function pass = test_or(pref)

if ( nargin < 1 )
    pref = chebfun.pref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Check the empty case.
f = chebfun(@sin);
pass(1) = isempty(f | chebfun()) && isempty(chebfun() | f);

% Check a few simple examples.
f = chebfun(@(x) sin(x), [-1 -0.5 0.5 1], pref);

g = chebfun(@(x) 0*x);
h = f | g;
ind = find(h.impulses == 0);
pass(2) = all(feval(h, x) == 1) && (numel(ind) == 1) && ...
    (abs(h.domain(ind)) < 10*vscale(h)*epslevel(h));

g = chebfun(@(x) exp(x));
h = f | g;
pass(3) = isequal(h.impulses, [1 1].') && all(feval(h, x) == 1);

f = chebfun(@(x) [sin(x) cos(x)], pref);
g = chebfun(@(x) [0*x exp(x)], pref);

h = f | g;
h_exact = @(x) [(sin(x) | 0*x) (0*x + 1)];
err = feval(h, x) - h_exact(x);
pass(4) = (norm(err(:), inf) == 0) && isequal(h.impulses, [1 1 ; 0 1 ; 1 1]);

h = f.' | g.';
h_exact = @(x) [(sin(x) | 0*x) (0*x + 1)].';
err = feval(h, x) - h_exact(x);
pass(5) = (norm(err(:), inf) == 0) && isequal(h.impulses, [1 1 ; 0 1 ; 1 1]);

% Check error conditions.
try
    f = chebfun(@(x) sin(x));
    g = chebfun(@(x) [sin(x) exp(x)]);
    h = f | g;
    pass(6) = false;
catch ME
    pass(6) = strcmp(ME.identifier, 'CHEBFUN:or:dims');
end

try
    f = chebfun(@(x) sin(x));
    g = chebfun(@(x) cos(x), [-2 7]);
    h = f | g;
    pass(7) = false;
catch ME
    pass(7) = strcmp(ME.identifier, 'CHEBFUN:or:doms');
end

end

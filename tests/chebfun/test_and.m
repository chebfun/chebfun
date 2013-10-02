% Test file for @chebfun/and.m.

function pass = test_and(pref)

if ( nargin < 1 )
    pref = chebfun.pref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Check the empty case.
f = chebfun(@sin);
pass(1) = isempty(f & chebfun()) && isempty(chebfun() & f);

% Check a few simple examples.
f = chebfun(@(x) sin(x), [-1 -0.5 0.5 1], pref);

g = chebfun(@(x) 0*x);
h = f & g;
pass(2) = all(h.impulses == 0) && all(feval(h, x) == 0);

g = chebfun(@(x) exp(x));
h = f & g;
ind = find(h.impulses == 0);
pass(3) = all(feval(h, x) == 1) && (numel(ind) == 1) && ...
    (abs(h.domain(ind)) < 10*vscale(h)*epslevel(h));

f = chebfun(@(x) [sin(x) cos(x)], pref);
g = chebfun(@(x) [exp(x) 0*x], pref);

h = f & g;
h_exact = @(x) [(sin(x) & exp(x)) 0*x];
err = feval(h, x) - h_exact(x);
pass(4) = (norm(err(:), inf) == 0) && isequal(h.impulses, [1 0 ; 0 0 ; 1 0]);

h = f.' & g.';
h_exact = @(x) [(sin(x) & exp(x)) 0*x].';
err = feval(h, x) - h_exact(x);
pass(5) = (norm(err(:), inf) == 0) && isequal(h.impulses, [1 0 ; 0 0 ; 1 0]);

% Check error conditions.
try
    f = chebfun(@(x) sin(x));
    g = chebfun(@(x) [sin(x) exp(x)]);
    h = f & g;
    pass(6) = false;
catch ME
    pass(6) = strcmp(ME.identifier, 'CHEBFUN:and:dims');
end

try
    f = chebfun(@(x) sin(x));
    g = chebfun(@(x) cos(x), [-2 7]);
    h = f & g;
    pass(7) = false;
catch ME
    pass(7) = strcmp(ME.identifier, 'CHEBFUN:and:doms');
end

end

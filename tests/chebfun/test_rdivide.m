% Test file for @chebfun/rdivide.m.

function pass = test_rdivide(pref)

if ( nargin < 1 )
    pref = chebfun.pref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Check the empty cases.
f = chebfun(@sin, pref);
g = chebfun();
pass(1) = isempty(f./g) && isempty(g./f);

% Check zero-numerator case.
pass(2) = iszero(0./f);

% Check a few simple examples.
f = chebfun(@(x) exp(x), [-1 -0.5 0 0.5 1], pref);
g = 1./f;
g_exact = @(x) exp(-x);
pass(3) = norm(feval(g, x) - g_exact(x), inf) < 10*vscale(g)*epslevel(g);

f = chebfun(@(x) exp(x), [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) exp(-x), [-1 1], pref);
h = f./g;
h_exact = @(x) exp(2*x);
pass(4) = norm(feval(h, x) - h_exact(x), inf) < 10*vscale(h)*epslevel(h);

f = chebfun(@(x) [exp(x) exp(-x)], [-1 -0.5 0 0.5 1], pref);
g = 1./f;
g_exact = @(x) [exp(-x) exp(x)];
err = feval(g, x) - g_exact(x);
pass(5) = norm(err(:), inf) < 10*vscale(g)*epslevel(g);

f = chebfun(@(x) [exp(x) exp(-x)], [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) [exp(-x) exp(x)], [-1 1], pref);
h = f./g;
h_exact = @(x) [exp(2*x) exp(-2*x)];
err = feval(h, x) - h_exact(x);
pass(6) = norm(err(:), inf) < 10*vscale(h)*epslevel(h);

ft = f.';
gt = g.';
h = ft./gt;
h_exact = @(x) [exp(2*x) exp(-2*x)].';
err = feval(h, x) - h_exact(x);
pass(7) = norm(err(:), inf) < 10*vscale(h)*epslevel(h);

% Check error conditions.
try
    h = f./0;
    pass(8) = false;
catch ME
    pass(8) = strcmp(ME.identifier, 'CHEBFUN:rdivide:DivisionByZero');
end

try
    f./gt;
    pass(9) = false;
catch ME
    pass(9) = strcmp(ME.identifier, 'CHEBFUN:rdivide:dim');
end

try
    f = chebfun(@(x) exp(x), [-1 1]);
    g = chebfun(@(x) exp(x), [0 2]);
    h = f./g;
    pass(10) = false;
catch ME
    pass(10) = strcmp(ME.identifier, 'CHEBFUN:rdivide:domain');
end

end

% Test file for @chebfun/fix.m

function pass = test_fix(pref)

if ( nargin < 1 )
    pref = chebpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Check empty case.
pass(1) = isempty(fix(chebfun()));

% Check "ordinary" cases.
f = chebfun(@(x) 2*sin(x), [-1 -0.5 0.5 1], pref);
g = fix(f);
g_exact = @(x) fix(2*sin(x));
pass(2) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*epslevel(g));

f = chebfun(@(x) exp(x), [-1 -0.5 0.5 1], pref);
g = fix(f);
g_exact = @(x) fix(exp(x));
pass(3) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*epslevel(g));

% Check complex-valued functions.
f2 = 1i*f;
g = fix(f2);
g_exact = @(x) fix(1i*exp(x));
pass(4) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*epslevel(g));

f3 = f + 1i*f;
g = fix(f3);
g_exact = @(x) fix(exp(x) + 1i*exp(x));
pass(5) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*epslevel(g));

% Check array-valued function.
f = chebfun(@(x) [2*sin(x) exp(x)], [-1 -0.5 0.5 1], pref);
g = fix(f);
g_exact = @(x) fix([2*sin(x) exp(x)]);
err = feval(g, x) - g_exact(x);
pass(6) = (norm(err(:), inf) <= 10*vscale(g)*epslevel(g));

% Check error conditions.
try
    f = chebfun(@(x) sin(x), [-1 -0.5 0 0.5 1], pref);
    f.impulses(3,1,2) = 1;
    g = fix(f);
    pass(7) = false;
catch ME
    pass(7) = strcmp(ME.identifier, 'CHEBFUN:fix:inf');
end

end
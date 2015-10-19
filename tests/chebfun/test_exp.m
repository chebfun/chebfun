% Test file for chebfun exp() and related functions.

function pass = test_exp(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

pref.splitting = 1;

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% List of functions to test.
expFunctions = {@exp, @expm1};

% Function with which we will be composing.
base_op = @(x) sign(x - 0.1).*abs(x + 0.2).*sin(3*x);
f = chebfun(base_op, [-1 -0.2 0.1 1], pref);

% Do the tests.
for (k = 1:1:numel(expFunctions))
    exp_op = expFunctions{k};
    g_exact = @(x) exp_op(base_op(x));
    g = exp_op(f, pref);
    err = feval(g, xr) - g_exact(xr);
    pass(k) = norm(err, inf) < 1e2*vscale(g)*eps;
    
end

end


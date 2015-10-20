function pass = test_besselj(pref)

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
xr = 2 * rand(100, 1) - 1;

% Arbitrary value to use for nu in the tests.
% [TODO]:  Use non-integer nu once fractional powers are implemented.
nu = 2;

f = chebfun(@(x) exp(x), [-1 -0.5 0 0.5 1], pref);

h = besselj(nu, f);
pass(1) = norm(feval(h, xr) - besselj(nu, exp(xr)), inf) < ...
    10*eps*vscale(h);

h2 = besselj(nu, f, 0);
pass(2) = normest(h - h2) < ...
    10*max(eps*vscale(h), eps*vscale(h2));

h = besselj(nu, f, 1);
pass(3) = norm(feval(h, xr) - besselj(nu, exp(xr), 1), inf) < ...
    10*eps*vscale(h);

% Test for array-valued chebfun.
f_op = @(x) [exp(-x) 1./(1 + 25*(x - 0.1).^2)];
f = chebfun(f_op, [-1 -0.5 0 0.5 1], pref);

h = besselj(nu, f);
err = feval(h, xr) - besselj(nu, f_op(xr));
pass(4) = norm(err(:), inf) < 1e2*eps*vscale(h);
    

h2 = besselj(nu, f, 0);
pass(5) = normest(h - h2) < ...
    10*max(eps*vscale(h), eps*vscale(h2));

h = besselj(nu, f, 1);
err = feval(h, xr) - besselj(nu, f_op(xr), 1);
pass(6) = norm(err(:), inf) < 1e2*eps*vscale(h);
    

%% Test for complex values.
pref.splitting = 1;
f_op = @complex_test_fn;
f = chebfun(f_op, [-1 0 0.5 1], pref);

h = besselj(nu, f, 0, pref);
pass(7) = norm(feval(h, xr) - besselj(nu, f_op(xr), 0), inf) < ...
    1e2*eps*vscale(h);

% Check for error on nu.
try
    f = chebfun(@(x) x, pref);
    h = besselj(1i + 3, f);
    pass(8) = false;
catch ME
    pass(8) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:besselj:nu');
end

end

function y = complex_test_fn(x)
    y = zeros(size(x));
    y(x <= 0) = exp(4*pi*1i*(x(x <= 0)));
    y(x > 0) = exp(x(x > 0));
end

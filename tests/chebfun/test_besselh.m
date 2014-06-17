% Test file for @chebfun/besselh.m.

function pass = test_besselh(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
xr = 2 * rand(100, 1) - 1;

% Enable splitting.

% Arbitrary value to use for nu in the tests.
nu = 1.663;

f = chebfun(@(x) exp(x), [-1 -0.5 0 0.5 1], pref);

h = besselh(nu, f);
pass(1) = norm(feval(h, xr) - besselh(nu, exp(xr)), inf) < ...
    10*epslevel(h)*vscale(h);

h2 = besselh(nu, 1, f);
pass(2) = normest(h - h2) < ...
    10*max(epslevel(h)*vscale(h), epslevel(h2)*vscale(h2));

h3 = besselh(nu, 1, f, 0);
pass(3) = normest(h - h3) < ...
    10*max(epslevel(h)*vscale(h), epslevel(h3)*vscale(h3));

h = besselh(nu, 2, f);
pass(4) = norm(feval(h, xr) - besselh(nu, 2, exp(xr)), inf) < ...
    10*epslevel(h)*vscale(h);

h = besselh(nu, 1, f, 1);
pass(5) = norm(feval(h, xr) - besselh(nu, 1, exp(xr), 1), inf) < ...
    10*epslevel(h)*vscale(h);

% Test for array-valued chebfun.
f_op = @(x) [exp(-x) 1./(1 + 25*(x - 0.1).^2)];
f = chebfun(f_op, [-1 -0.5 0 0.5 1], pref);

h = besselh(nu, f);
err = feval(h, xr) - besselh(nu, f_op(xr));
pass(6) = norm(err(:), inf) < 10*epslevel(h)*vscale(h);

h2 = besselh(nu, 1, f);
pass(7) = normest(h - h2) < ...
    10*max(epslevel(h)*vscale(h), epslevel(h2)*vscale(h2));

h3 = besselh(nu, 1, f, 0);
pass(8) = normest(h - h3) < ...
    10*max(epslevel(h)*vscale(h), epslevel(h3)*vscale(h3));

h = besselh(nu, 2, f);
err = feval(h, xr) - besselh(nu, 2, f_op(xr));
pass(9) = norm(err(:), inf) < 10*epslevel(h)*vscale(h);

h = besselh(nu, 1, f, 1);
err = feval(h, xr) - besselh(nu, 1, f_op(xr), 1);
pass(10) = norm(err(:), inf) < 10*epslevel(h)*vscale(h);

% Test for complex values.
pref.splitting = 1;
f_op = @scalar_test_fn;
f = chebfun(f_op, [-1 0 0.5 1], pref);

h = besselh(nu, 1, f, 0, pref);
pass(11) = norm(feval(h, xr) - besselh(nu, 1, f_op(xr), 0), inf) < ...
    10*epslevel(h)*vscale(h);

% Check for error on roots.
try
    f = chebfun(@(x) x, pref);
    h = besselh(nu, f);
    pass(12) = false;
catch ME
    pass(12) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:besselh:zero');
end

end

function y = scalar_test_fn(x)
    y = zeros(size(x));
    y(x <= 0) = exp(4*pi*1i*(x(x <= 0)));
    y(x > 0) = exp(x(x > 0));
end

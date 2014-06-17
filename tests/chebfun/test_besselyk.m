% Test file for @chebfun/bessely.m and @chebfun/besselk.m.

function pass = test_besselyk(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
xr = 2 * rand(100, 1) - 1;

% Arbitrary value to use for nu in the tests.
nu = 1.663;

for (n = 1:2)
    if (n == 1)
        testfn = @bessely;
        errid = 'CHEBFUN:CHEBFUN:bessely:zero';
    else
        testfn = @besselk;
        errid = 'CHEBFUN:CHEBFUN:besselk:zero';
    end

    f = chebfun(@(x) exp(x), [-1 -0.5 0 0.5 1], pref);

    h = testfn(nu, f);
    pass(n, 1) = norm(feval(h, xr) - testfn(nu, exp(xr)), inf) < ...
        10*epslevel(h)*vscale(h);

    h2 = testfn(nu, f, 0);
    pass(n, 2) = normest(h - h2) < ...
        10*max(epslevel(h)*vscale(h), epslevel(h2)*vscale(h2));

    h = testfn(nu, f, 1);
    pass(n, 3) = norm(feval(h, xr) - testfn(nu, exp(xr), 1), inf) < ...
        10*epslevel(h)*vscale(h);

    % Test for array-valued chebfun.
    f_op = @(x) [exp(-x) 1./(1 + 25*(x - 0.1).^2)];
    f = chebfun(f_op, [-1 -0.5 0 0.5 1], pref);

    h = testfn(nu, f);
    err = feval(h, xr) - testfn(nu, f_op(xr));
    pass(n, 4) = norm(err(:), inf) < 50*epslevel(h)*vscale(h);

    h2 = testfn(nu, f, 0);
    pass(n, 5) = normest(h - h2) < ...
        10*max(epslevel(h)*vscale(h), epslevel(h2)*vscale(h2));

    h = testfn(nu, f, 1);
    err = feval(h, xr) - testfn(nu, f_op(xr), 1);
    pass(n, 6) = norm(err(:), inf) < 10*epslevel(h)*vscale(h);

    % Test for complex values.
    pref.splitting = 1;
    f_op = @complex_test_fn;
    f = chebfun(f_op, [-1 0 0.5 1], pref);

    h = testfn(nu, f, 0, pref);
    pass(n, 7) = norm(feval(h, xr) - testfn(nu, f_op(xr), 0), inf) < ...
        10*epslevel(h)*vscale(h);

    % Check for error on roots.
    try
        f = chebfun(@(x) x, pref);
        h = testfn(nu, f);
        pass(n, 8) = false;
    catch ME
        pass(n, 8) = strcmp(ME.identifier, errid);
    end
end

end

function y = complex_test_fn(x)
    y = zeros(size(x));
    y(x <= 0) = exp(4*pi*1i*(x(x <= 0)));
    y(x > 0) = exp(x(x > 0));
end

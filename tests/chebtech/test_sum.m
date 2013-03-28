% Test file for chebtech/sum.

function pass = test_sum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.pref();
end

% Set a tolerance.
tol = 10*pref.chebtech.eps;

for ( n = 1:2 )
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    %%
    % Spot-check integrals for a couple of functions.
    f = testclass.make(@(x) exp(x) - 1, [], [], pref);
    pass(n, 1) = (abs(sum(f) - 0.350402387287603) < tol);

    f = testclass.make(@(x) 1./(1 + x.^2), [], [], pref);
    pass(n, 2) = (abs(sum(f) - pi/2) < tol);

    f = testclass.make(@(x) cos(1e4*x), [], [], pref);
    exact = -6.112287777765043e-05;
    pass(n, 3) = (abs(sum(f) - exact)/abs(exact) < 2e5*tol);

    z = exp(2*pi*1i/6);
    f = testclass.make(@(t) sinh(t*z), [], [], pref);
    pass(n, 4) = (abs(sum(f)) < tol);

    %%
    % Check a few basic properties.
    a = 2;
    b = -1i;
    f = testclass.make(@(x) x.*sin(x.^2) - 1, [], [], pref);
    df = diff(f);
    g = testclass.make(@(x) exp(-x.^2), [], [], pref);
    dg = diff(g);
    fg = f.*g;
    gdf = g.*df;
    fdg = f.*dg;

    % Linearity.
    pass(n, 5) = (abs(sum(a*f + b*g) - (a*sum(f) + b*sum(g))) < tol);

    % Integration-by-parts.
    pass(n, 6) = (abs(sum(fdg) - (feval(fg, 1) - feval(fg, -1) - sum(gdf))) ...
        < tol);

    % Fundamental Theorem of Calculus.
    pass(n, 7) = (abs(sum(df) - (feval(f, 1) - feval(f, -1))) < tol);
    pass(n, 8) = (abs(sum(dg) - (feval(g, 1) - feval(g, -1))) < tol);

    %%
    % Check operation for vectorized chebtech objects.
    f = testclass.make(@(x) [sin(x) x.^2 exp(1i*x)], [], [], pref);
    I = sum(f);
    I_exact = [0 2/3 2*sin(1)];
    pass(n, 9) = (max(abs(I - I_exact)) < tol);

    % Generate a few random points to use as test values.
    seedRNG(6178);
    x = 2 * rand(100, 1) - 1;

    % DIM option with vectorized input.
    g = sum(f, 2);
    h = @(x) sin(x) + x.^2 + exp(1i*x);
    pass(n, 10) = (norm(feval(g, x) - h(x), 'inf') < tol);

    % DIM option with non-vectorized input should leave everything alone.
    h = testclass.make(@(x) cos(x));
    sumh2 = sum(h, 2);
    pass(n, 11) = all((h.values == sumh2.values) & (h.coeffs == sumh2.coeffs));
end

end

% Test file for chebtech/sum.m

function pass = test_sum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.techPref();
end

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    %%
    % Spot-check integrals for a couple of functions.
    f = testclass.make(@(x) exp(x) - 1, [], pref);
    pass(n, 1) = (abs(sum(f) - 0.350402387287603) < 10*vscale(f)*eps);

    f = testclass.make(@(x) 1./(1 + x.^2), [], pref);
    pass(n, 2) = (abs(sum(f) - pi/2) < 10*vscale(f)*eps);

    f = testclass.make(@(x) cos(1e4*x), [], pref);
    exact = -6.112287777765043e-05;
    pass(n, 3) = (abs(sum(f) - exact)/abs(exact) < 1e6*vscale(f)*eps);
    
    
    z = exp(2*pi*1i/6);
    f = testclass.make(@(t) sinh(t*z), [], pref);
    pass(n, 4) = (abs(sum(f)) < 10*vscale(f)*eps);

    %%
    % Check a few basic properties.
    a = 2;
    b = -1i;
    f = testclass.make(@(x) x.*sin(x.^2) - 1, [], pref);
    df = diff(f);
    g = testclass.make(@(x) exp(-x.^2), [], pref);
    dg = diff(g);
    fg = f.*g;
    gdf = g.*df;
    fdg = f.*dg;

    tol_f = 10*vscale(f)*eps;
    tol_g = 10*vscale(f)*eps;
    tol_df = 10*vscale(df)*eps;
    tol_dg = 10*vscale(dg)*eps;
    tol_fg = 10*vscale(fg)*eps;
    tol_fdg = 10*vscale(fdg)*eps;
    tol_gdf = 10*vscale(gdf)*eps;

    % Linearity.
    pass(n, 5) = (abs(sum(a*f + b*g) - (a*sum(f) + b*sum(g))) < ...
        max(tol_f, tol_g));

    % Integration-by-parts.
    pass(n, 6) = (abs(sum(fdg) - (feval(fg, 1) - feval(fg, -1) - sum(gdf))) ...
        < max([tol_fdg ; tol_gdf ; tol_fg]));

    % Fundamental Theorem of Calculus.
    pass(n, 7) = (abs(sum(df) - (feval(f, 1) - feval(f, -1))) < ...
        max(tol_df, tol_f));
    pass(n, 8) = (abs(sum(dg) - (feval(g, 1) - feval(g, -1))) < ...
        max(tol_dg, tol_g));

    %%
    % Check operation for array-valued chebtech objects.
    f = testclass.make(@(x) [sin(x) x.^2 exp(1i*x)], [], pref);
    I = sum(f);
    I_exact = [0 2/3 2*sin(1)];
    pass(n, 9) = (max(abs(I - I_exact)) < 10*max(vscale(f)*eps));

    % Generate a few random points to use as test values.
    seedRNG(6178);
    x = 2 * rand(100, 1) - 1;

    % DIM option with array-valued input.
    g = sum(f, 2);
    h = @(x) sin(x) + x.^2 + exp(1i*x);
    pass(n, 10) = (norm(feval(g, x) - h(x), inf) < ...
        10*max(vscale(g)*eps));

    % DIM option with non-array-valued input should leave everything alone.
    h = testclass.make(@(x) cos(x));
    sumh2 = sum(h, 2);
    pass(n, 11) = all((h.coeffs == sumh2.coeffs));
end

end

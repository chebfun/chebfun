% Test file for chebtech/diff.m

function pass = test_diff(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    %%
    % Spot-check derivatives for a couple of functions.
    
    f = testclass.make(@(x) exp(x) - x, [], pref);
    df = diff(f);
    df_exact = @(x) exp(x) - 1;
    err = norm(df_exact(x) - feval(df, x), inf);
    tol = 100*df.vscale.*df.epslevel;
    pass(n, 1) = err < tol;
    
    f = testclass.make(@(x) atan(x), [], pref);
    df = diff(f);
    df_exact = @(x) 1./(1 + x.^2);
    err = norm(df_exact(x) - feval(df, x), inf);
    tol = 500*df.vscale.*df.epslevel;
    pass(n, 2) = err < tol;
    
    f = testclass.make(@(x) sin(x), [], pref);
    df = diff(f);
    df_exact = @(x) cos(x);
    err = norm(df_exact(x) - feval(df, x), inf);
    tol = 100*df.vscale.*df.epslevel;
    pass(n, 3) = err < tol;
    
    z = exp(2*pi*1i/3);
    f = testclass.make(@(t) airy(z*t), [], pref);
    df = diff(f);
    df_exact = @(t) z*airy(1, z*t);
    err = norm(df_exact(x) - feval(df, x), inf);
    tol = 100*df.vscale.*df.epslevel;
    pass(n, 4) = err < tol;
    
    %%
    % Verify that calling diff() gives the same answer as direct construction.
    
    f = testclass.make(@(x) 0.5*x - 0.0625*sin(8*x), [], pref);
    df = testclass.make(@(x) sin(4*x).^2, [], pref);
    err = diff(f) - df;
    pass(n, 5) = (norm(err.coeffs, inf) < 100*df.vscale.*df.epslevel);
    
    %%
    % Verify basic differentiation rules.
    
    f = testclass.make(@(x) x.*sin(x.^2) - 1, [], pref);
    df = diff(f);
    g = testclass.make(@(x) exp(-x.^2), [], pref);
    dg = diff(g);
    tol_f = 10*df.vscale.*df.epslevel;
    tol_g = 10*dg.vscale.*dg.epslevel;
    
    errfn = diff(f + g) - (df + dg);
    err = feval(errfn, x);
    pass(n, 6) = (norm(err, inf) < max(tol_f, tol_g));
    
    errfn = diff(f.*g) - (f.*dg + g.*df);
    err = feval(errfn, x);
    pass(n, 7) = (norm(err, inf) < length(f)*max(tol_f, tol_g));
    
    const = testclass.make(@(x) ones(size(x)), [], pref);
    dconst = diff(const);
    err = feval(dconst, x);
    pass(n, 8) = (norm(err, inf) == 0);
    
    %%
    % Check higher-order derivatives.  (NB:  We relax the tolerance by n + 1
    % factors of 10, where n is the number of derivatives taken.)
    
    f = testclass.make(@(x) x.*atan(x) - x - 0.5*log(1 + x.^2), [], pref);
    df2 = diff(f, 2);
    df2_exact = @(x) 1./(1 + x.^2);
    err = df2_exact(x) - feval(df2, x);
    pass(n, 9) = (norm(err, inf) < 300*df2.vscale.*df2.epslevel);
    
    f = testclass.make(@(x) sin(x), [], pref);
    df4 = diff(f, 4);
    df4_exact = @(x) sin(x);
    err = norm(df4_exact(x) - feval(df4, x), inf);
    tol = 900*df4.vscale.*df4.epslevel;
    pass(n, 10) = err < tol;

    f = testclass.make(@(x) x.^5 + 3*x.^3 - 2*x.^2 + 4, [], pref);
    df6 = diff(f, 6);
    df6_exact = @(x) zeros(size(x));
    err = df6_exact(x) - feval(df6, x);
    tol = 900*df4.vscale.*df4.epslevel;
    pass(n, 11) = (norm(err, inf) == 0);
    
    %%
    % Check operation for array-valued chebtech objects.
    
    f = testclass.make(@(x) [sin(x) x.^2 exp(1i*x)], [], pref);
    df = diff(f);
    df_exact = @(x) [cos(x) 2*x 1i*exp(1i*x)];
    err = feval(df, x) - df_exact(x);
    pass(n, 12) = (norm(err(:), inf) < 50*max(df.vscale.*df.epslevel));
    
    % DIM option.
    dim2df = diff(f, 1, 2);
    g = @(x) [(x.^2 - sin(x)) (exp(1i*x) - x.^2)];
    err = feval(dim2df, x) - g(x);
    pass(n, 13) = isequal(size(dim2df.vscale), [1 2]) && ...
        (norm(err(:), inf) < 10*max(dim2df.vscale.*dim2df.epslevel));

    dim2df2 = diff(f, 2, 2);
    g = @(x) exp(1i*x) - 2*x.^2 + sin(x);
    err = feval(dim2df2, x) - g(x);
    pass(n, 14) = isequal(size(dim2df2.vscale), [1 1]) && ...
        (norm(err(:), inf) < 10*max(dim2df2.vscale.*dim2df2.epslevel));

    % DIM option should return an empty chebtech for non-array-valued input.
    f = testclass.make(@(x) x.^3);
    dim2df = diff(f, 1, 2);
    pass(n, 15) = (isempty(dim2df.coeffs));
end

end

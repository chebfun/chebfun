% Test file for funcheb/diff.

function pass = test_diff(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb.pref;
end

% Set a tolerance.
tol = 1e3*pref.funcheb.eps;

% Generate a few random points to use as test values.
rngstate = rng();
rng(6178);
x = 2 * rand(100, 1) - 1;

for ( n = 1:2 )
    if ( n == 1 )
        testclass = funcheb1();
    else 
        testclass = funcheb2();
    end

    %%
    % Spot-check derivatives for a couple of functions.
    
    f = testclass.make(@(x) exp(x) - x, [], [], pref);
    df = diff(f);
    df_exact = @(x) exp(x) - 1;
    err = df_exact(x) - feval(df, x);
    pass(n, 1) = (norm(err, 'inf') < tol);
    
    f = testclass.make(@(x) atan(x), [], [], pref);
    df = diff(f);
    df_exact = @(x) 1./(1 + x.^2);
    err = df_exact(x) - feval(df, x);
    pass(n, 2) = (norm(err, 'inf') < tol);
    
    f = testclass.make(@(x) sin(x), [], [], pref);
    df = diff(f);
    df_exact = @(x) cos(x);
    err = df_exact(x) - feval(df, x);
    pass(n, 3) = (norm(err, 'inf') < tol);
    
    z = exp(2*pi*1i/3);
    f = testclass.make(@(t) airy(z*t), [], [], pref);
    df = diff(f);
    df_exact = @(t) z*airy(1, z*t);
    err = df_exact(x) - feval(df, x);
    pass(n, 4) = (norm(err, 'inf') < tol);
    
    %%
    % Verify that calling diff() gives the same answer as direct construction.
    
    f = testclass.make(@(x) 0.5*x - 0.0625*sin(8*x), [], [], pref);
    df = testclass.make(@(x) sin(4*x).^2, [], [], pref);
    err = diff(f) - df;
    pass(n, 5) = (norm(err.values, 'inf') < tol);
    
    %%
    % Verify basic differentiation rules.
    
    f = testclass.make(@(x) x.*sin(x.^2) - 1, [], [], pref);
    df = diff(f);
    g = testclass.make(@(x) exp(-x.^2), [], [], pref);
    dg = diff(g);
    
    errfn = diff(f + g) - (df + dg);
    err = feval(errfn, x);
    pass(n, 6) = (norm(err, 'inf') < tol);
    
    errfn = diff(f.*g) - (f.*dg + g.*df);
    err = feval(errfn, x);
    pass(n, 7) = (norm(err, 'inf') < length(f)*tol);
    
    const = testclass.make(@(x) ones(size(x)), [], [], pref);
    dconst = diff(const);
    err = feval(dconst, x);
    pass(n, 8) = (norm(err, 'inf') < tol);
    
    %%
    % Check higher-order derivatives.  (NB:  We relax the tolerance by n + 1
    % factors of 10, where n is the number of derivatives taken.)
    
    f = testclass.make(@(x) x.*atan(x) - x - 0.5*log(1 + x.^2), [], [], pref);
    df2 = diff(f, 2);
    df2_exact = @(x) 1./(1 + x.^2);
    err = df2_exact(x) - feval(df2, x);
    pass(n, 9) = (norm(err, 'inf') < 1e3*tol);
    
    f = testclass.make(@(x) sin(x), [], [], pref);
    df4 = diff(f, 4);
    df4_exact = @(x) sin(x);
    err = df4_exact(x) - feval(df4, x);
    pass(n, 10) = (norm(err, 'inf') < 1e5*tol);
    
    %%
    % Check operation for vectorized funcheb objects.
    
    f = testclass.make(@(x) [sin(x) x.^2 exp(1i*x)], [], [], pref);
    df_exact = @(x) [cos(x) 2*x 1i*exp(1i*x)];
    err = feval(diff(f), x) - df_exact(x);
    pass(n, 11) = (norm(err(:), 'inf') < tol);
    
    % DIM option.
    dim2df = diff(f, 1, 2);
    g = @(x) [(x.^2 - sin(x)) (exp(1i*x) - x.^2)];
    err = feval(dim2df, x) - g(x);
    pass(n, 12) = (norm(err(:), 'inf') < tol);
    
    dim2df2 = diff(f, 2, 2);
    g = @(x) exp(1i*x) - 2*x.^2 + sin(x);
    err = feval(dim2df2, x) - g(x);
    pass(n, 13) = (norm(err(:), 'inf') < tol);
    
    % DIM option should return an empty funcheb for non-vectorized input.
    f = testclass.make(@(x) x.^3);
    dim2df = diff(f, 1, 2);
    pass(n, 14) = (isempty(dim2df.values) && isempty(dim2df.coeffs));
end

%%
% Restore the RNG state.
rng(rngstate);

end


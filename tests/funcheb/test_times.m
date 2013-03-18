% Test file for funcheb/times.

function pass = test_times(pref)

% Get preferences.
if (nargin < 1)
    pref = funcheb.pref;
end

% Generate a few random points to use as test values.
rngstate = rng();
rng(6178);
x = 2 * rand(100, 1) - 1;

% Random numbers to use as arbitrary multiplicative constants.
alpha = -0.194758928283640 + 0.075474485412665i;
beta = -0.526634844879922 - 0.685484380523668i;

for ( n = 1:2 )
    if ( n == 1 )
        testclass = funcheb1();
    else 
        testclass = funcheb2();
    end

    %%
    % Check operation in the face of empty arguments.
    
    f = testclass.make();
    g = testclass.make(@(x) x, [], [], pref);
    pass(n, 1) = (isempty(f .* f) && isempty(f .* g) && isempty(g .* f));
    
    %%
    % Check multiplication by scalars.
    
    f_op = @(x) sin(x);
    f = testclass.make(f_op, [], [], pref);
    pass(n, 2:3) = test_mult_function_by_scalar(f, f_op, alpha, x);
    
    f_op = @(x) [sin(x) cos(x)];
    f = testclass.make(f_op, [], [], pref);
    pass(n, 4:5) = test_mult_function_by_scalar(f, f_op, alpha, x);
    
    % Can't multiply by matrix of scalars with more than one row.
    try 
        disp(f .* [1 2 ; 3 4]);
        pass(n, 6) = false;
    catch ME
        pass(n, 6) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:times:dim');
    end
    
    % This should fail with a dimension mismatch error.
    try
        disp(f .* [1 2 3]);
        pass(n, 7) = false;
    catch ME
        pass(n, 7) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:times:dim');
    end
    
    %%
    % Check multiplication by constant functions.
    
    f_op = @(x) sin(x);
    f = testclass.make(f_op, [], [], pref);
    g_op = @(x) alpha*ones(size(x));
    g = testclass.make(g_op, [], [], pref);
    pass(n, 8) = test_mult_function_by_function(f, f_op, g, g_op, x, false);
    
    % This should fail with a dimension mismatch error from funcheb.mtimes().
    f_op = @(x) [sin(x) cos(x)];
    f = testclass.make(f_op, [], [], pref);
    g_op = @(x) repmat([alpha, beta], size(x, 1), 1);
    g = testclass.make(g_op, [], [], pref);
    pass(n, 9) = test_mult_function_by_function(f, f_op, g, g_op, x, false);
    
    %%
    % Spot-check multiplication of two funcheb objects for a few test 
    % functions.
    
    f_op = @(x) ones(size(x));
    f = testclass.make(f_op, [], [], pref);
    pass(n, 10) = test_mult_function_by_function(f, f_op, f, f_op, x, false);
    
    f_op = @(x) exp(x) - 1;
    f = testclass.make(f_op, [], [], pref);
    
    g_op = @(x) 1./(1 + x.^2);
    g = testclass.make(g_op, [], [], pref);
    pass(n, 11) = test_mult_function_by_function(f, f_op, g, g_op, x, false);
    
    g_op = @(x) cos(1e4*x);
    g = testclass.make(g_op, [], [], pref);
    pass(n, 12) = test_mult_function_by_function(f, f_op, g, g_op, x, false);
    
    g_op = @(t) sinh(t*exp(2*pi*1i/6));
    g = testclass.make(g_op, [], [], pref);
    pass(n, 13) = test_mult_function_by_function(f, f_op, g, g_op, x, false);
    
    %%
    % Check operation for vectorized funcheb objects.
    
    f = testclass.make(@(x) [sin(x) cos(x) exp(x)], [], [], pref);
    g = testclass.make(@(x) tanh(x), [], [], pref);
    h1 = f .* g;
    h2 = g .* f;
    pass(n, 14) = isequal(h1, h2);
    h_exact = @(x) [tanh(x).*sin(x) tanh(x).*cos(x) tanh(x).*exp(x)];
    err = feval(h1, x) - h_exact(x);
    pass(n, 15) = max(abs(err(:))) < 10*h1.epslevel;
    
    g = testclass.make(@(x) [sinh(x) cosh(x) tanh(x)], [], [], pref);
    h = f .* g;
    h_exact = @(x) [sinh(x).*sin(x) cosh(x).*cos(x) tanh(x).*exp(x)];
    err = feval(h, x) - h_exact(x);
    pass(n, 16) = max(abs(err(:))) < 10*h.epslevel;
    
    % This should fail with a dimension mismatch error.
    try
        g = testclass.make(@(x) [sinh(x) cosh(x)], [], [], pref);
        disp(f .* g);
        pass(n, 17) = false;
    catch ME
        pass(n, 17) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:times:dim2');
    end
    
    %%
    % Check specially handled cases, including some in which an adjustment for
    % positivity is performed.
    
    f_op = @(t) sinh(t*exp(2*pi*1i/6));
    f = testclass.make(f_op, [], [], pref);
    pass(n, 18) = test_mult_function_by_function(f, f_op, f, f_op, x, false);
    
    g_op = @(t) conj(sinh(t*exp(2*pi*1i/6)));
    g = conj(f);
    pass(n, 19:20) = test_mult_function_by_function(f, f_op, g, g_op, x, true);
    
    f_op = @(x) exp(x) - 1;
    f = testclass.make(f_op, [], [], pref);
    pass(n, 21:22) = test_mult_function_by_function(f, f_op, f, f_op, x, true);
    
    %%
    % Check that multiplication and direct construction give similar results.
    
    tol = 10*eps;
    g_op = @(x) 1./(1 + x.^2);
    g = testclass.make(g_op, [], [], pref);
    h1 = f .* g;
    h2 = testclass.make(@(x) f_op(x) .* g_op(x), [], [], pref);
    h2 = prolong(h2, length(h1));
    pass(n, 23) = norm(h1.values - h2.values, 'inf') < tol;
end

%%
% Restore the RNG state.

rng(rngstate);

end

% Test the multiplication of a FUNCHEB F, specified by F_OP, by a scalar ALPHA
% using a grid of points X in [-1  1] for testing samples.
function result = test_mult_function_by_scalar(f, f_op, alpha, x)
    g1 = f .* alpha;
    g2 = alpha .* f;
    result(1) = isequal(g1, g2);
    g_exact = @(x) f_op(x) .* alpha;
    result(2) = norm(feval(g1, x) - g_exact(x), 'inf') < 10*g1.epslevel;
end

% Test the multiplication of two FUNCHEB objects F and G, specified by F_OP and
% G_OP, using a grid of points X in [-1  1] for testing samples.  If CHECKPOS is
% TRUE, an additional check is performed to ensure that the values of the result
% are all nonnegative; otherwise, this check is skipped.
function result = test_mult_function_by_function(f, f_op, g, g_op, x, checkpos)
    h = f .* g;
    h_exact = @(x) f_op(x) .* g_op(x);
    result(1) = norm(feval(h, x) - h_exact(x), 'inf') < 10*h.epslevel;
    if ( checkpos )
        result(2) = all(h.values >= 0);
    end
end

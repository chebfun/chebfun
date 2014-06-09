% Test file for chebtech/plus.m

function pass = test_plus(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% A random number to use as an arbitrary additive constant.
alpha = -0.194758928283640 + 0.075474485412665i;

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    %%
    % Check operation in the face of empty arguments.
    
    f = testclass.make();
    g = testclass.make(@(x) x, [], pref);
    pass(n, 1) = (isempty(f + f) && isempty(f + g) && isempty(g + f));
    
    %%
    % Check addition with scalars.
    
    f_op = @(x) sin(x);
    f = testclass.make(f_op, [], pref);
    pass(n, 2:3) = test_add_function_to_scalar(f, f_op, alpha, x);
    
    %%
    % Check addition of two chebtech objects.
    
    f_op = @(x) zeros(size(x));
    f = testclass.make(f_op, [], pref);
    pass(n, 4:5) = test_add_function_to_function(f, f_op, f, f_op, x);
    
    f_op = @(x) exp(x) - 1;
    f = testclass.make(f_op, [], pref);
    
    g_op = @(x) 1./(1 + x.^2);
    g = testclass.make(g_op, [], pref);
    pass(n, 6:7) = test_add_function_to_function(f, f_op, g, g_op, x);
    
    g_op = @(x) cos(1e4*x);
    g = testclass.make(g_op, [], pref);
    pass(n, 8:9) = test_add_function_to_function(f, f_op, g, g_op, x);
    
    g_op = @(t) sinh(t*exp(2*pi*1i/6));
    g = testclass.make(g_op, [], pref);
    pass(n, 10:11) = test_add_function_to_function(f, f_op, g, g_op, x);
    
    %%
    % Check operation for array-valued chebtech objects.
    
    f_op = @(x) [zeros(size(x)) zeros(size(x)) zeros(size(x))];
    f = testclass.make(f_op, [], pref);
    pass(n, 12:13) = test_add_function_to_function(f, f_op, f, f_op, x);
    
    f_op = @(x) [sin(x) cos(x) exp(x)];
    f = testclass.make(f_op, [], pref);
    pass(n, 14:15) = test_add_function_to_scalar(f, f_op, alpha, x);
    
    g_op = @(x) [cosh(x) airy(1i*x) sinh(x)];
    g = testclass.make(g_op, [], pref);
    pass(n, 16:17) = test_add_function_to_function(f, f_op, g, g_op, x);
    
    % This should fail with a dimension mismatch error.
    g_op = @(x) sin(x);
    g = testclass.make(g_op, [], pref);
    try
        h = f + g; %#ok<NASGU>
        pass(n, 18) = false;
    catch ME
        pass(n, 18) = strcmp(ME.message, 'Matrix dimensions must agree.');
    end
    
    %%
    % Check that direct construction and PLUS give comparable results.
    
    tol = 10*eps;
    f = testclass.make(@(x) x, [], pref);
    g = testclass.make(@(x) cos(x) - 1, [], pref);
    h1 = f + g;
    h2 = testclass.make(@(x) x + cos(x) - 1, [], pref);
    
    pass(n, 19) = norm(h1.coeffs - h2.coeffs, inf) < tol;

    %%
    % Check that adding a CHEBTECH to an unhappy CHEBTECH gives an unhappy
    % result.  

    f = testclass.make(@(x) cos(x+1));    % Happy
    g = testclass.make(@(x) sqrt(x+1));   % Unhappy
    h = f + g;  % Add unhappy to happy.
    pass(n, 20) = (~g.ishappy) && (~h.ishappy);
    h = g + f;  % Add happy to unhappy.
    pass(n, 21) = (~g.ishappy) && (~h.ishappy);

    %%
    % Test addition of array-valued scalar to array-valued chebtech.

    f = testclass.make(@(x) [sin(x) cos(x) exp(x)]);
    g = f + [1 2 3];
    g_exact = @(x) [(1 + sin(x)) (2 + cos(x)) (3 + exp(x))];
    err = feval(g, x) - g_exact(x);
    pass(n, 22) = norm(err(:), inf) < 10*max(g.vscale.*g.epslevel);

    %%
    % Test scalar expansion in CHEBTECH argument.

    f = testclass.make(@(x) sin(x));
    g = f + [1 2 3];
    g_exact = @(x) [(1 + sin(x)) (2 + sin(x)) (3 + sin(x))];
    err = feval(g, x) - g_exact(x);
    pass(n, 23) = isequal(size(g.coeffs, 2), 3) && norm(err(:), inf) < ...
        10*max(g.vscale.*g.epslevel);
end

end

% Test the addition of a CHEBTECH F, specified by F_OP, to a scalar ALPHA using
% a grid of points X in [-1  1] for testing samples.
function result = test_add_function_to_scalar(f, f_op, alpha, x)
    g1 = f + alpha;
    g2 = alpha + f;
    result(1) = isequal(g1, g2);
    g_exact = @(x) f_op(x) + alpha;
    result(2) = norm(feval(g1, x) - g_exact(x), inf) <= ...
        10*max(g1.vscale.*g1.epslevel);
end

% Test the addition of two CHEBTECH objects F and G, specified by F_OP and
% G_OP, using a grid of points X in [-1  1] for testing samples.
function result = test_add_function_to_function(f, f_op, g, g_op, x)
    h1 = f + g;
    h2 = g + f;
    result(1) = isequal(h1, h2);
    h_exact = @(x) f_op(x) + g_op(x);
    norm(feval(h1, x) - h_exact(x), inf);
    result(2) = norm(feval(h1, x) - h_exact(x), inf) <= ...
        10*max(h1.vscale.*h1.epslevel);
end

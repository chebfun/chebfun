% Test file for chebtech/rdivide.m

function pass = test_rdivide(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Random numbers to use as arbitrary constants.
alpha = -0.194758928283640 + 0.075474485412665i;
beta = -0.526634844879922 - 0.685484380523668i;

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    %%
    % Check division by single scalars.
    
    f_op = @(x) sin(x);
    f = testclass.make(f_op, [], pref);
    pass(n, 1) = test_div_function_by_scalar(f, f_op, alpha, x);
    
    g = f ./ 0;
    pass(n, 2) = isnan(g);
    
    f_op = @(x) [sin(x) cos(x)];
    f = testclass.make(f_op, [], pref);
    pass(n, 3) = test_div_function_by_scalar(f, f_op, alpha, x);
    
    g = f ./ 0;
    pass(n, 4) = isnan(g);
    
    %%
    % Check division by a row matrix of scalars in the case of an array-valued
    % chebtech object.
    
    g = f ./ [alpha beta];
    g_exact = @(x) [sin(x)./alpha cos(x)./beta];
    pass(n, 5) = norm(feval(g, x) - g_exact(x), inf) < 10*eps;
    
    g = f ./ [alpha 0];
    pass(n, 6) = isnan(g) && ~any(isnan(g.coeffs(:, 1))) ...
                       && all(isnan(g.coeffs(:, 2)));
    
    %%
    % Check division of a scalar by a chebtech object.
    
    f_op = @(x) exp(x);
    f = testclass.make(@(x) exp(x), [], pref);
    pass(n, 7) = test_div_scalar_by_function(alpha, f, f_op, x);
    
    %%
    % Check division of two chebtech objects.
    
    g_op = @(x) exp(x);
    g = testclass.make(g_op, [], pref);
    
    f_op = @(x) exp(x) - 1;
    f = testclass.make(f_op, [], pref);
    pass(n, 8) = test_div_function_by_function(f, f_op, g, g_op, x);
    
    f_op = @(x) 1./(1 + x.^2);
    f = testclass.make(f_op, [], pref);
    pass(n, 9) = test_div_function_by_function(f, f_op, g, g_op, x);
    
    f_op = @(x) cos(1e4*x);
    f = testclass.make(f_op, [], pref);
    pass(n, 10) = test_div_function_by_function(f, f_op, g, g_op, x);
    
    f_op = @(t) sinh(t*exp(2*pi*1i/6));
    f = testclass.make(f_op, [], pref);
    pass(n, 11) = test_div_function_by_function(f, f_op, g, g_op, x);
    
    %%
    % Check proper behavior under error conditions.
    
    % Can't divide by a scalar matrix with multiple rows.
    try
        f = testclass.make(@(x) sin(x), [], pref);
        disp(f ./ [1 ; 2]);
        pass(n, 12) = false;
    catch ME
        pass(n, 12) = strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:rdivide:size');
    end
    
    % Can't divide by a scalar row matrix if the column counts don't match.
    try
        f = testclass.make(@(x) [sin(x) cos(x)], [], pref);
        disp(f ./ [1 2 3]);
        pass(n, 13) = false;
    catch ME
        pass(n, 13) = strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:rdivide:size');
    end
    
    %%
    % Check that direct construction and RDIVIDE give comparable results.
    
    tol = 100*eps;
    
    f = testclass.make(@(x) sin(x), [], pref);
    h1 = f ./ alpha;
    h2 = testclass.make(@(x) sin(x) ./ alpha, [], pref);
    pass(n, 14) = norm(h1.coeffs - h2.coeffs, inf) < tol;
    
    g = testclass.make(@(x) exp(x), [], pref);
    h1 = f ./ g;
    h2 = testclass.make(@(x) sin(x) ./ exp(x), [], pref);
    maxlengthtmp = max(length(h1),length(h2));
    h1 = prolong(h1,maxlengthtmp);
    h2 = prolong(h2,maxlengthtmp);
    pass(n, 15) = norm(h1.coeffs - h2.coeffs, inf) < tol;
end

end

% Test the division of a CHEBTECH F, specified by F_OP, by a scalar ALPHA using
% a grid of points X in [-1  1] for testing samples.
function result = test_div_function_by_scalar(f, f_op, alpha, x)
    g = f ./ alpha;
    g_exact = @(x) f_op(x) ./ alpha;
    result = norm(feval(g, x) - g_exact(x), inf) < 10*max(vscale(g)*eps);
end

% Test the division of a scalar ALPHA by a CHEBTECH, specified by F_OP, using
% a grid of points X in [-1  1] for testing samples.
function result = test_div_scalar_by_function(alpha, f, f_op, x)
    g = alpha ./ f;
    g_exact = @(x) alpha ./ f_op(x);
    result = norm(feval(g, x) - g_exact(x), inf) < 10*max(vscale(g)*eps);
end

% Test the division of two CHEBTECH objects F and G, specified by F_OP and
% G_OP, using a grid of points X in [-1  1] for testing samples.
function result = test_div_function_by_function(f, f_op, g, g_op, x)
    h = f ./ g;
    h_exact = @(x) f_op(x) ./ g_op(x);
    norm(feval(h, x) - h_exact(x), inf);
    result = norm(feval(h, x) - h_exact(x), inf) < 1e4*max(vscale(h)*eps);
        
        % (1e2 is enough except solely for test 10, which requires bigger)
end

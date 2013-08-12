% Test file for fun/rdivide.m

function pass = test_rdivide(pref)

% Get preferences.
if ( nargin < 1 )
    pref = fun.pref;
end

% Generate a few random points to use as test values.
seedRNG(6178);

% Random numbers to use as arbitrary constants.
alpha = -0.194758928283640 + 0.075474485412665i;
beta = -0.526634844879922 - 0.685484380523668i;

pass = zeros(1, 13); % Pre-allocate pass matrix
for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
        
        % Set the domain
        dom = [-2 7];
        x = diff(dom) * rand(100, 1) + dom(1);
    else 
        testclass = unbndfun();
    end

    %%
    % Check division by single scalars.
    
    f_op = @(x) sin(x);
    f = testclass.make(f_op, dom, [], [], pref);
    pass(n, 1) = test_div_function_by_scalar(f, f_op, alpha, x);
    
    g = f ./ 0;
    pass(n, 2) = isnan(g);
    
    f_op = @(x) [sin(x) cos(x)];
    f = testclass.make(f_op, dom, [], [], pref);
    pass(n, 3) = test_div_function_by_scalar(f, f_op, alpha, x);
    
    g = f ./ 0;
    pass(n, 4) = isnan(g);
    
    %%
    % Check division by a row matrix of scalars in the case of an array-valued
    % fun object.
    
    g = f ./ [alpha beta];
    g_exact = @(x) [sin(x)./alpha cos(x)./beta];
    pass(n, 5) = norm(feval(g, x) - g_exact(x), inf) < ...
        10*max(get(g, 'vscale'))*get(g, 'epslevel');
    
    g = f ./ [alpha 0];
    isn = isnan(feval(g, x));
    pass(n, 6) = isnan(g) && ~any(any(isn(:,1))) ...
                       && all(isn(:,2));
    
    %%
    % Check division of a scalar by a fun object.
    
    f_op = @(x) exp(x);
    f = testclass.make(@(x) exp(x), dom, [], [], pref);
    pass(n, 7) = test_div_scalar_by_function(alpha, f, f_op, x);
    
    %%
    % Check division of two fun objects.
    
    g_op = @(x) exp(x);
    g = testclass.make(g_op, dom, [], [], pref);
    
    f_op = @(x) exp(x) - 1;
    f = testclass.make(f_op, dom, [], [], pref);
    pass(n, 8) = test_div_function_by_function(f, f_op, g, g_op, x);
    
    f_op = @(x) 1./(1 + x.^2);
    f = testclass.make(f_op, dom, [], [], pref);
    pass(n, 9) = test_div_function_by_function(f, f_op, g, g_op, x);
    
    f_op = @(x) cos(1e4*x);
    f = testclass.make(f_op, dom, [], [], pref);
    pass(n, 10) = test_div_function_by_function(f, f_op, g, g_op, x);
    
    f_op = @(t) sinh(t*exp(2*pi*1i/6));
    f = testclass.make(f_op, dom, [], [], pref);
    pass(n, 11) = test_div_function_by_function(f, f_op, g, g_op, x);
    
    %%
    % Check that direct construction and RDIVIDE give comparable results.
    
    f = testclass.make(@(x) sin(x), dom, [], [], pref);
    h1 = f ./ alpha;
    h2 = testclass.make(@(x) sin(x) ./ alpha, dom, [], [], pref);
    pass(n, 12) = norm(feval(h1, x) - feval(h2, x), inf) < ...
        10*get(h2, 'vscale')*get(h2, 'epslevel');
    
    g = testclass.make(@(x) exp(x), dom, [], [], pref);
    h1 = f ./ g;
    h2 = testclass.make(@(x) sin(x) ./ exp(x), dom, [], [], pref);
    pass(n, 13) = norm(feval(h1, x) - feval(h2, x), inf) < ...
        5e3*get(h2, 'vscale')*get(h2, 'epslevel');
end

end

% Test the division of a FUN F, specified by F_OP, by a scalar ALPHA using
% a grid of points X in [-1  1] for testing samples.
function result = test_div_function_by_scalar(f, f_op, alpha, x)
    g = f ./ alpha;
    g_exact = @(x) f_op(x) ./ alpha;
    result = norm(feval(g, x) - g_exact(x), inf) < ...
        10*max(get(g, 'vscale'))*get(g, 'epslevel');
end

% Test the division of a scalar ALPHA by a FUN, specified by F_OP, using
% a grid of points X in [-1  1] for testing samples.
function result = test_div_scalar_by_function(alpha, f, f_op, x)
    g = alpha ./ f;
    g_exact = @(x) alpha ./ f_op(x);
    result = norm(feval(g, x) - g_exact(x), inf) < ...
        10*max(get(g, 'vscale'))*get(g, 'epslevel');
end

% Test the division of two FUN objects F and G, specified by F_OP and
% G_OP, using a grid of points X in [-1  1] for testing samples.
function result = test_div_function_by_function(f, f_op, g, g_op, x)
    h = f ./ g;
    h_exact = @(x) f_op(x) ./ g_op(x);
    norm(feval(h, x) - h_exact(x), inf);
    result = norm(feval(h, x) - h_exact(x), inf) < ...
        50*max(get(h, 'vscale'))*get(h, 'epslevel');
end

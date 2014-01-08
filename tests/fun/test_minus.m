% Test file for fun/minus.m

function pass = test_minus(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);

% A random number to use as an arbitrary additive constant.
alpha = randn() + 1i*randn();

for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
        dom = [-2 7];
        x = diff(dom) * rand(100, 1) + dom(1);
    else 
        testclass = unbndfun();
    end

    %%
    % Check operation in the face of empty arguments.
    
    f = testclass.make();
    g = testclass.make(@(x) x, dom, [], [], pref);
    pass(n, 1) = (isempty(f - f) && isempty(f - g) && isempty(g - f));
    
    %%
    % Check subtraction with scalars.
    
    f_op = @(x) sin(x);
    f = testclass.make(f_op, dom, [], [], pref);
    pass(n, 2:3) = test_sub_function_and_scalar(f, f_op, alpha, x);
    
    %%
    % Check subtraction of two bndfun objects.
    
    f_op = @(x) zeros(size(x));
    f = testclass.make(f_op, dom, [], [], pref);
    pass(n, 4:5) = test_sub_function_and_function(f, f_op, f, f_op, x);
    
    f_op = @(x) exp(x) - 1;
    f = testclass.make(f_op, dom, [], [], pref);
    
    g_op = @(x) 1./(1 + x.^2);
    g = testclass.make(g_op, dom, [], [], pref);
    pass(n, 6:7) = test_sub_function_and_function(f, f_op, g, g_op, x);
    
    g_op = @(x) cos(1e4*x);
    g = testclass.make(g_op, dom, [], [], pref);
    pass(n, 8:9) = test_sub_function_and_function(f, f_op, g, g_op, x);
    
    g_op = @(t) sinh(t*exp(2*pi*1i/6));
    g = testclass.make(g_op, dom, [], [], pref);
    pass(n, 10:11) = test_sub_function_and_function(f, f_op, g, g_op, x);
    
    %%
    % Check operation for array-valued bndfun objects.
    
    f_op = @(x) [zeros(size(x)) zeros(size(x)) zeros(size(x))];
    f = testclass.make(f_op, dom, [], [], pref);
    pass(n, 12:13) = test_sub_function_and_function(f, f_op, f, f_op, x);
    
    f_op = @(x) [sin(x) cos(x) exp(x)];
    f = testclass.make(f_op, dom, [], [], pref);
    pass(n, 14:15) = test_sub_function_and_scalar(f, f_op, alpha, x);
    
    g_op = @(x) [cosh(x) airy(1i*x) sinh(x)];
    g = testclass.make(g_op, dom, [], [], pref);
    pass(n, 16:17) = test_sub_function_and_function(f, f_op, g, g_op, x);
    
    % This should fail with a dimension mismatch error.
    g_op = @(x) sin(x);
    g = testclass.make(g_op, dom, [], [], pref);
    try
        h = f - g; %#ok<NASGU>
        pass(n, 18) = false;
    catch ME
        pass(n, 18) = strcmp(ME.message, 'Matrix dimensions must agree.');
    end
    
    %%
    % Check that direct construction and MINUS give comparable results.
    
    tol = 10*eps;
    f = testclass.make(@(x) x, dom, [], [], pref);
    g = testclass.make(@(x) cos(x) - 1, dom, [], [], pref);
    h1 = f - g;
    h2 = testclass.make(@(x) x - (cos(x) - 1), dom, [], [], pref);
    pass(n, 19) = normest(h1-h2) < tol;

    %%
    % Check that subtracting a BNDFUN and an unhappy BNDFUN gives an
    % unhappy result.  

    f = testclass.make(@(x) cos(x+1), dom);    % Happy
    g = testclass.make(@(x) sqrt(x+1), dom);   % Unhappy
    h = f - g;  % Subtract unhappy from happy.
    pass(n, 20) = (~get(g, 'ishappy')) && (~get(h, 'ishappy'));
    h = g - f;  % Subtract happy from unhappy.
    pass(n, 21) = (~get(g, 'ishappy')) && (~get(h, 'ishappy'));
end

end

% Test the subtraction of a BNDFUN F, specified by F_OP, to and from a scalar
% ALPHA using a grid of points X in [a  b] for testing samples.
function result = test_sub_function_and_scalar(f, f_op, alpha, x)
    g1 = f - alpha;
    g2 = alpha - f;
    result(1) = isequal(g1, -g2);
    g_exact = @(x) f_op(x) - alpha;
    result(2) = norm(feval(g1, x) - g_exact(x), inf) < ...
        10*max(get(f,'vscale').*get(g1, 'epslevel'));
end

% Test the subraction of two BNDFUN objects F and G, specified by F_OP and
% G_OP, using a grid of points X in [-1  1] for testing samples.
function result = test_sub_function_and_function(f, f_op, g, g_op, x)
    h1 = f - g;
    h2 = g - f;
    result(1) = isequal(h1, -h2);
    h_exact = @(x) f_op(x) - g_op(x);
    norm(feval(h1, x) - h_exact(x), inf);
    result(2) = norm(feval(h1, x) - h_exact(x), inf) <= 10*max(get(h1,'vscale').*get(h1,'epslevel'));
end

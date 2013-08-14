% Test file for fun/mtimes.m

function pass = test_mtimes(pref)

% Get preferences.
if ( nargin < 1 )
    pref = fun.pref;
end

% Create seed for random number generator
seedRNG(6178);

% A random number to use as an arbitrary scalar multiplier.
alpha = randn() + 1i*randn();

pass = zeros(1, 12); % Pre-allocate pass matrix
for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
        dom = [-2 7];
        
        % Generate a few random points to use as test values.
        x = diff(dom) * rand(1000, 1) + dom(1);
        
    else
        testclass = unbndfun();
    end

    %%
    % Check operation in the face of empty arguments.
    
    f = testclass.make(@(x) sin(x), dom, [], [], pref);
    g = testclass.make();
    pass(n, 1) = isempty(f*[]) && isempty([]*f) && isempty(2*g) && isempty(g*2);
    
    %%
    % Check operation for scalar fun objects.
    
    f = testclass.make(@(x) sin(x), dom, [], [], pref);
    g1 = alpha*f;
    g2 = f*alpha;
    pass(n, 2) = isequal(g1, g2);
    g_exact = @(x) alpha*sin(x);
    pass(n, 3) = norm(feval(g1, x) - g_exact(x), inf) < get(g1, 'vscale')*get(g1, 'epslevel');
    
    g = 0*f;
    pass(n, 4) = all(feval(g, x) == 0);
    
    %%
    % Check operation for array-valued fun objects.
    
    f = testclass.make(@(x) [sin(x) cos(x) exp(x)], dom, [], [], pref);
    g1 = alpha*f;
    g2 = f*alpha;
    pass(n, 5) = isequal(g1, g2);
    g_exact = @(x) alpha*[sin(x) cos(x) exp(x)];
    err = abs(feval(g1, x) - g_exact(x));
    pass(n, 6) = max(err(:)) < max(get(g1, 'vscale'))*get(g1, 'epslevel');
    
    g = 0*f;
    pass(n, 7) = all(all(feval(g, x) == zeros(numel(x), 3)));
    
    A = randn(3, 3);
    g = f*A;
    g_exact = @(x) [sin(x) cos(x) exp(x)]*A;
    err = abs(feval(g, x) - g_exact(x));
    pass(n, 8) = max(err(:)) < 2*max(get(g, 'vscale'))*get(g, 'epslevel');
    
    %%
    % Verify error handling and corner cases.
    
    % Multiply non-scalar double and fun.
    try
        f = testclass.make(@(x) exp(x), dom);
        disp([1 2 3]*f)
    catch ME
        pass(n, 9) = strcmp(ME.identifier, 'CHEBFUN:FUN:mtimes:size') ...
            && strcmp(ME.message, 'Inner matrix dimensions must agree.');
    end
    
    % Multiply fun and non-scalar double with mismatching dimensions.
    try
        f = testclass.make(@(x) [sin(x) cos(x)], dom);
        disp(f*[1 ; 2 ; 3]);
        pass(n, 10) = false;
    catch ME
        pass(n, 10) = strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:mtimes:size2') ...
            && strcmp(ME.message, 'Inner matrix dimensions must agree.');
    end
    
    % Using * for multiplication of two fun objects.
    try
        g = testclass.make(@(x) x, dom);
        disp(f*g);
        pass(n, 11) = false;
    catch ME
        pass(n, 11) = strcmp(ME.message, 'Use .* to multiply FUN objects.');
    end
    
    % Using * to multiply a fun and something else.
    try
        disp(f*uint8(128));
        pass(n, 12) = false;
    catch ME
        pass(n, 12) = strcmp(ME.message, ...
            'mtimes does not know how to multiply a FUN and a uint8.');
    end
end

end

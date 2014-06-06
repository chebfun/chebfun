% Test file for chebtech/mrdivide.m

function pass = test_mrdivide(pref)

% Get preferences.
if (nargin < 1)
    pref = chebtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Random number to use as a scalar constant.
alpha = -0.194758928283640 + 0.075474485412665i;

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    %%
    % Check division of a chebtech object by a numeric array.
    
    f_op = @(x) [sin(x) cos(x)];
    f = testclass.make(f_op, [], pref);
    pass(n, 1) = isnan(f / 0);
    
    g = f / alpha;
    g_exact = @(x) f_op(x) ./ alpha;
    pass(n, 2) = norm(feval(g, x) - g_exact(x), inf) < ...
        10*max(g.vscale.*g.epslevel);
    
    % A "least-squares" case where the solution is obvious.
    I = eye(2);
    g = f / I;
    err = g*I - f;
    pass(n, 3) = max(max(abs(feval(err, x)))) < 10*max(g.vscale.*g.epslevel);
    
    % A less trivial least-squares case for which we still know the answer.
    A = [1 1];
    g = f / A;
    g_exact = @(x) (sin(x) + cos(x))/2;
    pass(n, 4) = norm(feval(g, x) - g_exact(x), inf) < ...
        10*max(g.vscale.*g.epslevel);
    
    %%
    % Check division of a numeric array by a chebtech object.
    
    f = testclass.make(@(x) sin(x));
    g = alpha / f;
    pass(n, 5) = abs(innerProduct(f, g) - alpha) < 10*max(g.vscale.*g.epslevel);
    
    f = testclass.make(@(x) [sin(2*pi*x) cos(2*pi*x)]);
    g = [1 1]/f;
    g_exact = @(x) (sin(2*pi*x) + cos(2*pi*x));
    pass(n, 6) = norm(feval(g, x) - g_exact(x), inf) < ...
        10*max(g.vscale.*g.epslevel);
    
    %%
    % Check error conditions.
    
    % Catch dimension mismatch errors.
    try
        g = f / [1 2 3];
        pass(n, 7) = false;
    catch ME
        pass(n, 7) = strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:mrdivide:size');
    end
    
    % Can't do f/g if both f and g are chebtech objects.
    try
        f = testclass.make(@(x) sin(x));
        g = testclass.make(@(x) cos(x));
        h = f / g;
        pass(n, 8) = false;
    catch ME
        pass(n, 8) = strcmp(ME.identifier, ...
            'CHEBFUN:CHEBTECH:mrdivide:chebtechDivChebtech');
    end
    
    % Can't call mrdivide on a chebtech and a non-chebtech or non-double 
    % object.
    try
        f = testclass.make(@(x) sin(x));
        g = f / true;
        pass(n, 9) = false;
    catch ME
        pass(n, 9) = strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:mrdivide:badArg');
    end
end

end

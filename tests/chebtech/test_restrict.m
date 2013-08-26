% Test file for chebtech/restrict.m

function pass = test_restrict(pref)

% Get preferences.
if (nargin < 1)
    pref = chebtech.pref;
end

pass = zeros(2, 12); % Pre-allocate pass matrix
for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else
        testclass = chebtech2();
    end

    %%
    % Check behavior for empty inputs.
    f = testclass.make();
    f = restrict(f, [-0.5 0.5]);
    pass(n, 1) = isempty(f);

    %%
    % Check behvaior for non-subinterval inputs.
    f = testclass.make(@(x) sin(x), [], [], pref);
    g = restrict(f, [-1 1]);
    pass(n, 2) = isequal(f, g);

    try
        g = restrict(f, [-1, 3]); %#ok<NASGU>
        pass(n, 3) = 0;
    catch ME
        pass(n, 3) = strcmp(ME.identifier, 'CHEBTECH:restrict:badinterval');
    end

    try
        g = restrict(f, [-2, 1]); %#ok<NASGU>
        pass(n, 4) = 0;
    catch ME
        pass(n, 4) = strcmp(ME.identifier, 'CHEBTECH:restrict:badinterval');
    end

    try
        g = restrict(f, [-1 -0.25 0.3 0.1 1]); %#ok<NASGU>
        pass(n, 5) = 0;
    catch ME
        pass(n, 5) = strcmp(ME.identifier, 'CHEBTECH:restrict:badinterval');
    end

    %%
    % Spot-check a few functions
    pass(n, 6) = test_spotcheck_restrict(testclass, ...
        @(x) exp(x) - 1, [-0.2 0.1], pref);
    pass(n, 7) = test_spotcheck_restrict(testclass, ...
        @(x) 1./(1 + x.^2), [-0.7 0.9], pref);
    pass(n, 8) = test_spotcheck_restrict(testclass, ...
        @(x) cos(1e3*x), [0.1 0.5], pref);
    pass(n, 9) = test_spotcheck_restrict(testclass, ...
        @(t) sinh(t*exp(2*pi*1i/6)), [-0.4 1], pref);

    %%
    % Check multiple subinterval restriction.
    f = testclass.make(@(x) sin(x) + sin(x.^2), [], [], pref);
    g = restrict(f, [-0.7 0.3 0.8]);
    h1 = restrict(f, [-0.7 0.3]);
    h2 = restrict(f, [0.3 0.8]);
    x = linspace(-1, 1, 100).';
    err1 = feval(g{1} - h1, x);
    err2 = feval(g{2} - h2, x);
    pass(n, 10) = all(err1(:) == 0) && all(err2(:) == 0);

    %%
    % Check operation for array-valued functions.
    pass(n, 11) = test_spotcheck_restrict(testclass, ...
        @(x) [sin(x) cos(x) exp(x)], [-1 -0.7], pref);

    f = testclass.make(@(x) [sin(x) cos(x)], [], [], pref);
    g = restrict(f, [-0.6 0.1 1]);
    h1 = restrict(f, [-0.6 0.1]);
    h2 = restrict(f, [0.1 1]);
    x = linspace(-1, 1, 100).';
    err1 = feval(g{1} - h1, x);
    err2 = feval(g{2} - h2, x);
    pass(n, 12) = all(err1(:) == 0) && all(err2(:) == 0);
end

end

% Spot-check restriction of a given function to a given subinterval.
function result = test_spotcheck_restrict(testclass, fun_op, subint,pref)
    % Perform restriction.
    f = testclass.make(fun_op, [], [], pref);
    g = restrict(f, subint);

    % Construct mapping from restricted subinterval to [-1, 1].
    a = subint(1);
    b = subint(2);
    map = @(t) (2/(b - a))*(t - a) - 1;

    % Sample on a grid of 100 points and check for accuracy.
    x = linspace(a, b, 100).';
    y_exact = fun_op(x);
    y_approx = feval(g, map(x));

    result = norm(y_exact - y_approx, Inf) < 10*max(g.vscale.*g.epslevel);
end

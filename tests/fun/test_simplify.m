% Test file for fun/simplify.m

function pass = test_simplify(pref)

% Get preferences:
if ( nargin < 1 )
    pref = fun.pref;
end

% Set the domain
dom = [-1 1];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

pass = zeros(1, 11); % Pre-allocate pass matrix
for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
        % Tolerance for passing to simplify:
        simptol = 1e-6;
    else
        testclass = unbndfun();
    end
   
    %%
    % Test pathological inputs.
    
    % Empty FUN objects should be left alone.
    f = testclass.make();
    g = simplify(f);
    pass(n, 1) = isequal(f, g);

    % Unhappy FUN objects should be left alone.
    f = testclass.make(@(x) sqrt(x), dom);
    g = simplify(f);
    pass(n, 2) = ~f.onefun.ishappy && isequal(f, g);

    %%
    % Test for a scalar-valued function:

    f = testclass.make(@(x) sin(100*(x + 0.1)), dom);
    g = simplify(f, simptol);
    pass(n, 3) = length(g) < length(f);
    pass(n, 4) = norm(feval(f, x) - feval(g, x), inf) < 10*get(g,'epslevel')*get(g,'vscale');

    %%
    % Lengths of simplifications should be invariant under scaling:

    f1 = 1e-8*f;
    g1 = simplify(f1, simptol);
    pass(n, 5) = length(g1) == length(g);

    f2 = 1e8*f;
    g2 = simplify(f2, simptol);
    pass(n, 6) = length(g2) == length(g);

    %%
    % Test for an array-valued function:

    f = testclass.make(@(x) [sin(100*(x + 0.1)) cos(100*(x + 0.1)) exp(x)], dom);
    g = simplify(f, 1e-6);
    pass(n, 7) = length(g) < length(f);
    pass(n, 8) = all(norm(feval(f, x) - feval(g, x), inf) < ...
        10*get(g,'epslevel')*get(g,'vscale'));

    %%
    % Try a contrived example which will leave a zero FUN:

    f = testclass.make(@(x) sin(100*(x + 0.1)), dom);
    g = simplify(f, 1e20);
    pass(n, 9) = iszero(g);

    %%
    % Try an example that zeros only interior coefficients, not the tail:

    f = testclass.make(@(x) (1) + (1e-10*x) + (4*x.^3 - 3*x), dom);
    g = simplify(f, 1e-6);
    pass(n, 10) = length(g) == length(f);
    pass(n, 11) = norm(feval(f, x) - feval(g, x), inf) < 10*get(g,'epslevel')*get(g,'vscale');
end

end


% Test file for chebtech/simplify.m

function pass = test_simplify(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Tolerance for passing to simplify:
simptol = 1e-6;

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else
        testclass = chebtech2();
    end

    %%
    % Test pathological inputs.

    % Empty CHEBTECH objects should be left alone.
    f = testclass.make();
    g = simplify(f);
    pass(n, 1) = isequal(f, g);

    % Unhappy CHEBTECH objects should be left alone.
    f = testclass.make(@(x) sqrt(x));
    g = simplify(f);
    pass(n, 2) = ~f.ishappy && isequal(f, g);

    %%
    % Test for a scalar-valued function:

    f = testclass.make(@(x) sin(100*(x + 0.1)));
    g = simplify(f, simptol);
    pass(n, 3) = abs(g.coeffs(end)) ~= 0;
    pass(n, 4) = length(g) < length(f);
    pass(n, 5) = norm(feval(f, x) - feval(g, x), inf) < 1e2*simptol*vscale(f);

    %%
    % Lengths of simplifications should be invariant under scaling:

    f1 = 1e-8*f;
    g1 = simplify(f1, simptol);
    pass(n, 6) = all(abs(g1.coeffs) >= eps*vscale(g1));
    pass(n, 7) = length(g1) == length(g);

    f2 = 1e8*f;
    g2 = simplify(f2, simptol);
    pass(n, 8) = all(abs(g2.coeffs) >= eps*vscale(g2));
    pass(n, 9) = length(g2) == length(g);

    %%
    % Test for an array-valued function:

    f = testclass.make(@(x) [sin(100*(x + 0.1)) cos(100*(x + 0.1)) exp(x)]);
    g = simplify(f, simptol);
    pass(n, 10) = any(abs(g.coeffs(1, :)) ~= 0);
    pass(n, 11) = length(g) < length(f);
    pass(n, 12) = all(norm(feval(f, x) - feval(g, x), inf) < ...
        10*max(simptol.*vscale(f)));

    %%
    % Try a contrived example which will leave a length 1 CHEBTECH:

    f = testclass.make(@(x) sin(100*(x + 0.1)));
    g = simplify(f, 1e20);
    pass(n, 13) = length(g) == 1;

    % Check that a long identically-zero CHEBTECH simplifies correctly:
    f = testclass.make(@(x) 0*x, [], struct('fixedLength', 8));
    g = simplify(f);
    pass(n, 14) = iszero(g) && (length(g) == 1);

end

end

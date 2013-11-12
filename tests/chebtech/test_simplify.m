% Test file for chebtech/simplify.m

function pass = test_simplify(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else
        testclass = chebtech2();
    end

    % Tolerance for passing to simplify:
    simptol = 1e-6;

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
    pass(n, 3) = all((abs(g.coeffs) > simptol*g.vscale) | (g.coeffs == 0));
    pass(n, 4) = abs(g.coeffs(end)) ~= 0;
    pass(n, 5) = length(g) < length(f);
    pass(n, 6) = norm(feval(f, x) - feval(g, x), inf) < 10*g.epslevel*g.vscale;

    %%
    % Lengths of simplifications should be invariant under scaling:

    f1 = 1e-8*f;
    g1 = simplify(f1, simptol);
    pass(n, 7) = all(abs(g1.coeffs) >= g1.epslevel*g1.vscale);
    pass(n, 8) = length(g1) == length(g);

    f2 = 1e8*f;
    g2 = simplify(f2, simptol);
    pass(n, 9) = all(abs(g2.coeffs) >= g1.epslevel*g2.vscale);
    pass(n, 10) = length(g2) == length(g);

    %%
    % Test for an array-valued function:

    f = testclass.make(@(x) [sin(100*(x + 0.1)) cos(100*(x + 0.1)) exp(x)]);
    g = simplify(f, 1e-6);
    pass(n, 11) = all(all((abs(g.coeffs) > ...
        repmat(simptol*g.vscale, length(g), 1)) | (g.coeffs == 0)));
    pass(n, 12) = any(abs(g.coeffs(1, :)) ~= 0);
    pass(n, 13) = length(g) < length(f);
    pass(n, 14) = all(norm(feval(f, x) - feval(g, x), inf) < ...
        10*g.epslevel*g.vscale);

    %%
    % Try a contrived example which will leave a zero CHEBTECH:

    f = testclass.make(@(x) sin(100*(x + 0.1)));
    g = simplify(f, 1e20);
    pass(n, 15) = iszero(g);

    %%
    % Try an example that zeros only interior coefficients, not the tail:

    f = testclass.make(@(x) (1) + (1e-10*x) + (4*x.^3 - 3*x));
    g = simplify(f, 1e-6);
    pass(n, 16) = all((abs(g.coeffs) > simptol*g.vscale) | (g.coeffs == 0));
    pass(n, 17) = abs(g.coeffs(end)) ~= 0;
    pass(n, 18) = length(g) == length(f);
    pass(n, 19) = norm(feval(f, x) - feval(g, x), inf) < 10*g.epslevel*g.vscale;
end

end

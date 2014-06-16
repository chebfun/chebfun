% Test file for chebtech/prolong.m

function pass = test_prolong(pref)

if ( nargin < 1 )
    pref = chebtech.techPref();
end

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    F = @sin;
    f = testclass.make(@(x) F(x), [], pref);

    k = 11;
    g = prolong(f, k);
    x = testclass.chebpts(k);
    values = g.coeffs2vals(g.coeffs);
    pass(n, 1) = norm(values - F(x), inf) < 10*g.vscale.*g.epslevel;

    g = prolong(f, 1);
    pass(n, 2) = size(g,1) == 1 && norm(g.coeffs, inf) < ...
        10*g.vscale.*g.epslevel;

    g = prolong(f, 2);
    exact_values = sin(testclass.chebpts(2));
    values = g.coeffs2vals(g.coeffs);
    pass(n, 3) = size(g,1) == 2 && norm(values - exact_values, inf) < ...
        10*g.vscale.*g.epslevel;

    F = @sin;
    f = testclass.make(@(x) [F(x), -F(x)], [], pref);
    k = 11;
    g = prolong(f, k);
    x = testclass.chebpts(k);
    values = g.coeffs2vals(g.coeffs);
    pass(n, 4) = norm(values - [F(x), -F(x)], inf) < ...
        10*max(g.vscale.*g.epslevel);

    g = prolong(f, 1);
    pass(n, 5) = size(g,1) == 1 && norm(g.coeffs, inf) < ...
        10*max(g.vscale.*g.epslevel);

    g = prolong(f, 2);
    y = testclass.chebpts(2);
    exact_values = [sin(y) -sin(y)];
    values = g.coeffs2vals(g.coeffs);
    pass(n, 6) = size(g,1) == 2 && ...
        norm(values - exact_values, inf) < 10*max(g.vscale.*g.epslevel);

    g = prolong(f, length(f));
    fvalues = f.coeffs2vals(f.coeffs);
    gvalues = g.coeffs2vals(g.coeffs);
    pass(n, 7) = all(fvalues(:) == gvalues(:));

    k = 32;
    g = prolong(f, k);
    x = testclass.chebpts(k);
    values = g.coeffs2vals(g.coeffs);
    pass(n, 8) = norm(values - [F(x), -F(x)], inf) < ...
        10*max(g.vscale.*g.epslevel);

    k = 100;
    g = prolong(f, k);
    x = testclass.chebpts(k);
    values = g.coeffs2vals(g.coeffs);
    pass(n, 9) = norm(values - [F(x), -F(x)], inf) < ...
        10*max(g.vscale.*g.epslevel);

    F = @(x) sin(1000*x);
    f = testclass.make(@(x) [F(x), -F(x)], [], pref);
    k = 32;
    g = prolong(f, k);
    x = testclass.chebpts(k);
    values = g.coeffs2vals(g.coeffs);
    pass(n, 10) = norm(values - [F(x), -F(x)], inf) < ...
        10*max(g.vscale.*g.epslevel);

    k = 100;
    g = prolong(f, k);
    x = testclass.chebpts(k);
    values = g.coeffs2vals(g.coeffs);
    pass(n, 11) = norm(values - [F(x), -F(x)], inf) < ...
        10*max(g.vscale.*g.epslevel);

    g = prolong(f, 1);
    pass(n, 12) = size(g, 1) == 1 && norm(g.coeffs, inf) < ...
        10*max(g.vscale.*g.epslevel);

    g = prolong(f, 2);
    y = testclass.chebpts(2);
    exact_values = [sin(1000*y) -sin(1000*y)];
    values = g.coeffs2vals(g.coeffs);
    pass(n, 13) = size(g, 1) == 2 && ...
        norm(values - exact_values, inf) < 10*max(g.vscale.*g.epslevel);

    F = @(x) cos(1000*x);
    f = testclass.make(@(x) [F(x), -F(x)], [], pref);

    g = prolong(f, 1);
    values = g.coeffs2vals(g.coeffs);
    pass(n, 14) = size(g, 1) == 1 && norm(values - [1, -1], inf) < ...
        10*max(g.vscale.*g.epslevel);

    g = prolong(f, 2);
    y = testclass.chebpts(2);
    exact_values = [cos(1000*y) -cos(1000*y)];
    values = g.coeffs2vals(g.coeffs);
    pass(n, 15) = length(g) == 2 && ...
        norm(values - exact_values, inf) < 10*max(g.vscale.*g.epslevel);

    v = [1 2 3];
    f = testclass.make(v, [], pref);
    g = prolong(f, 5);
    values = g.coeffs2vals(g.coeffs);
    pass(n, 16) = norm(values - repmat([1 2 3], 5, 1), inf) < ...
       10*max(g.vscale.*g.epslevel);
end

end

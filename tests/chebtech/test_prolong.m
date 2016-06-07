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

    g = prolong(f, 1);
    pass(n, 1) = size(g,1) == 1 && norm(g.coeffs, inf) < 10*eps;

    F = @sin;
    f = testclass.make(@(x) [F(x), -F(x)], [], pref);
    g = prolong(f, 1);
    pass(n, 2) = size(g,1) == 1 && norm(g.coeffs, inf) < 10*eps;

    g = prolong(f, length(f));
    fvalues = f.coeffs2vals(f.coeffs);
    gvalues = g.coeffs2vals(g.coeffs);
    pass(n, 3) = all(fvalues(:) == gvalues(:));

    k = 32;
    g = prolong(f, k);
    x = testclass.chebpts(k);
    values = g.coeffs2vals(g.coeffs);
    pass(n, 4) = norm(values - [F(x), -F(x)], inf) < 10*max(vscale(g)*eps);

    k = 100;
    g = prolong(f, k);
    x = testclass.chebpts(k);
    values = g.coeffs2vals(g.coeffs);
    pass(n, 5) = norm(values - [F(x), -F(x)], inf) < 10*max(vscale(g)*eps);

    F = @(x) sin(1000*x);
    f = testclass.make(@(x) [F(x), -F(x)], [], pref);
    g = prolong(f, 1);
    pass(n, 6) = size(g, 1) == 1 && norm(g.coeffs, inf) < 1e2*eps;

    v = [1 2 3];
    f = testclass.make(v, [], pref);
    g = prolong(f, 5);
    values = g.coeffs2vals(g.coeffs);
    pass(n, 7) = norm(values - repmat([1 2 3], 5, 1), inf) < ...
       10*max(vscale(g)*eps);
end

end

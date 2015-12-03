% Test file for chebtech/alias.m

function pass = test_alias(pref)

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
    g = f;
    g.coeffs = testclass.alias(f.coeffs, k);
    x = testclass.chebpts(k);
    values = testclass.coeffs2vals(g.coeffs);
    pass(n, 1) = norm(values - F(x), inf) < 10*vscale(g)*eps;

    g.coeffs = testclass.alias(f.coeffs, 2);
    exact_values = sin(testclass.chebpts(2));
    values = testclass.coeffs2vals(g.coeffs);
    pass(n, 2) = size(g,1) == 2 && norm(values - exact_values, inf) < ...
        1e1*vscale(g)*eps;

    F = @sin;
    f = testclass.make(@(x) [F(x), -F(x)], [], pref);
    k = 11;
    g = f;
    g.coeffs = testclass.alias(f.coeffs, k);
    x = testclass.chebpts(k);
    values = testclass.coeffs2vals(g.coeffs);
    pass(n, 3) = norm(values - [F(x), -F(x)], inf) < 10*max(vscale(g)*eps);

    g = f;
    g.coeffs = testclass.alias(f.coeffs, 2);
    y = testclass.chebpts(2);
    exact_values = [sin(y) -sin(y)];
    values = testclass.coeffs2vals(g.coeffs);
    pass(n, 4) = size(g,1) == 2 && ...
        norm(values - exact_values, inf) < 10*max(vscale(g)*eps);

    F = @(x) sin(1000*x);
    f = testclass.make(@(x) [F(x), -F(x)], [], pref);
    k = 32;
    g = f;
    g.coeffs = testclass.alias(f.coeffs, k);
    x = testclass.chebpts(k);
    values = testclass.coeffs2vals(g.coeffs);
    pass(n, 5) = norm(values - [F(x), -F(x)], inf) < 1e3*max(vscale(g)*eps);
    
    k = 100;
    g = f;
    g.coeffs = testclass.alias(f.coeffs, k);
    x = testclass.chebpts(k);
    values = testclass.coeffs2vals(g.coeffs);
    pass(n, 6) = norm(values - [F(x), -F(x)], inf) < 1e4*max(vscale(g)*eps);
    
    g = f;
    g.coeffs = testclass.alias(f.coeffs, 2);
    y = testclass.chebpts(2);
    exact_values = [sin(1000*y) -sin(1000*y)];
    values = testclass.coeffs2vals(g.coeffs);
    pass(n, 7) = size(g, 1) == 2 && ...
        norm(values - exact_values, inf) < 1e3*eps;
    

    F = @(x) cos(1000*x);
    f = testclass.make(@(x) [F(x), -F(x)], [], pref);

    g = f;
    g.coeffs = testclass.alias(f.coeffs, 1);
    values = testclass.coeffs2vals(g.coeffs);
    pass(n, 8) = size(g, 1) == 1 && norm(values - [1, -1], inf) < ...
        1e2*max(vscale(g)*eps);

    g = f;
    g.coeffs = testclass.alias(f.coeffs, 2);
    y = testclass.chebpts(2);
    exact_values = [cos(1000*y) -cos(1000*y)];
    values = testclass.coeffs2vals(g.coeffs);
    pass(n, 9) = length(g) == 2 && ...
        norm(values - exact_values, inf) < 1e3*max(vscale(g)*eps);
    
end

end

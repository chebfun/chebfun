function pass = test_prolong(pref)

if ( nargin < 1 )
    pref = funcheb.pref;
end
tol = 10*pref.funcheb.eps;

for ( n = 1:2 )
    if ( n == 1 )
        testclass = funcheb1();
    else 
        testclass = funcheb2();
    end

    F = @sin;
    f = testclass.make(@(x) F(x), [], [], pref);

    k = 11;
    g = prolong(f, k);
    x = testclass.chebpts(k);
    pass(n, 1) = norm(g.values - F(x), inf) < tol;

    g = prolong(f,1);
    pass(n, 2) = size(g,1) == 1 && norm(g.values, inf) < tol;

    g = prolong(f,2);
    exact_values = sin(testclass.chebpts(2));
    pass(n, 3) = size(g,1) == 2 && norm(g.values - exact_values, inf) < tol;

    F = @sin;
    f = testclass.make(@(x) [F(x), -F(x)], [], [], pref);
    k = 11;
    g = prolong(f, k);
    x = testclass.chebpts(k);
    pass(n, 4) = norm(g.values - [F(x), -F(x)], inf) < tol;

    g = prolong(f,1);
    pass(n, 5) = size(g,1) == 1 && norm(g.values, inf) < tol;

    g = prolong(f,2);
    y = testclass.chebpts(2);
    exact_values = [sin(y) -sin(y)];
    pass(n, 6) = size(g,1) == 2 && ...
        norm(g.values - exact_values, inf) < tol;

    g = prolong(f, length(f));
    pass(n, 7) = all(f.values(:) == g.values(:));

    k = 32;
    g = prolong(f, k);
    x = testclass.chebpts(k);
    pass(n, 8) = norm(g.values - [F(x), -F(x)], inf) < tol;

    k = 100;
    g = prolong(f, k);
    x = testclass.chebpts(k);
    pass(n, 9) = norm(g.values - [F(x), -F(x)], inf) < tol;

    F = @(x) sin(1000*x);
    f = testclass.make(@(x) [F(x), -F(x)], [], [], pref);
    k = 32;
    g = prolong(f, k);
    x = testclass.chebpts(k);
    pass(n, 10) = norm(g.values - [F(x), -F(x)], inf) < length(f)*tol;

    k = 100;
    g = prolong(f, k);
    x = testclass.chebpts(k);
    pass(n, 11) = norm(g.values - [F(x), -F(x)], inf) < length(f)*tol;

    g = prolong(f,1);
    pass(n, 12) = size(g, 1) == 1 && norm(g.values, inf) < length(f)*tol;

    g = prolong(f,2);
    y = testclass.chebpts(2);
    exact_values = [sin(1000*y) -sin(1000*y)];
    pass(n, 13) = size(g, 1) == 2 && ...
        norm(g.values - exact_values, inf) < length(f)*tol;

    F = @(x) cos(1000*x);
    f = testclass.make(@(x) [F(x), -F(x)], [], [], pref);

    g = prolong(f,1);
    pass(n, 14) = size(g, 1) == 1 && norm(g.values - [1, -1], inf) < tol;

    g = prolong(f,2);
    y = testclass.chebpts(2);
    exact_values = [cos(1000*y) -cos(1000*y)];
    pass(n, 15) = length(g) == 2 && ...
        norm(g.values - exact_values, inf) < length(f)*tol;

    v = [1 2 3];
    f = testclass.make(v, [], [], pref);
    g = prolong(f, 5);
    pass(n, 16) = norm(g.values - repmat([1 2 3], 5, 1), inf) < tol;
end

end

function pass = test_flipud(pref)

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

    f = testclass.make(@(x) sin(x+.5), 0, pref);
    g = testclass.make(@(x) sin(-x+.5), 0, pref);
    h = flipud(f);
    pass(n, 1) = norm(g.values - h.values, inf) < tol;
    
    f = testclass.make(@(x) [sin(x+.5), exp(x)], 0, pref);
    g = testclass.make(@(x) [sin(-x+.5), exp(-x)], 0, pref);
    h = flipud(f);
    pass(n, 2) = norm(g.values - h.values, inf) < tol;
    
    f = testclass.make(@(x) sin(1i*x+.5), 0, pref);
    g = testclass.make(@(x) sin(-1i*x+.5), 0, pref);
    h = flipud(f);
    pass(n, 3) = norm(g.values - h.values, inf) < tol;
    
    f = testclass.make(@(x) [sin(x+.5), exp(1i*x)], 0, pref);
    g = testclass.make(@(x) [sin(-x+.5), exp(-1i*x)], 0, pref);
    h = flipud(f);
    pass(n, 4) = norm(g.values - h.values, inf) < tol;
end

end

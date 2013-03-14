function pass = test_isfinite(pref)

if ( nargin < 1 )
    pref = funcheb.pref;
end
p = pref;

for ( n = 1:2 )
    if ( n == 1 )
        testclass = funcheb1();
    else 
        testclass = funcheb2();
    end

    % Test a scalar-valued function:
    p.funcheb.n = 11;
    y = testclass.chebpts(p.funcheb.n); % Force singularity to fall on grid.
    f = testclass.make(@(x) 1./(x - y(4)), 0, p);
    pass(n, 1) = ~isfinite(f);
    
    % Test a vector-valued function:
    p.funcheb.n = 11;
    y = testclass.chebpts(p.funcheb.n); % Force singularity to fall on grid.
    f = testclass.make(@(x) [1./(x - y(4)), x], 0, p);
    pass(n, 2) = ~isfinite(f);
    
    % Test a finite scalar-valued function:
    p = pref;
    f = testclass.make(@(x) x, 0, p);
    pass(n, 3) = isfinite(f);
    
    % Test a finite vector-valued function:
    f = testclass.make(@(x) [x, x], 0, p);
    pass(n, 4) = isfinite(f);
end

end

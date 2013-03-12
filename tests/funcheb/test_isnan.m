function pass = test_isnan(pref)

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
    f = testclass.make(@(x) x, p);
    pass(n, 1) = ~isnan(f);
    
    % Test a vector-valued function:
    f = testclass.make(@(x) [x, x.^2], p);
    pass(n, 2) = ~isnan(f);
    
    % Test a NaN scalar-valued function:
    try
        f = testclass.make(@(x) x + NaN);
        pass(n, 3) = isnan(f);
    catch ME
        pass(n, 3) = strcmpi(ME.message, 'Too many NaNs to handle.');
    end
    
    % Test a NaN vector-valued function:
    try
        f = testclass.make(@(x) [x, x + NaN]);
        pass(n, 4) = isnan(f);
    catch ME
        pass(n, 4) = strcmpi(ME.message, 'Too many NaNs to handle.');
    end
    
    % Test a NaN vector-valued function:
    % [TODO]:  This test fails for funcheb1.
    try
        f = testclass.make(@(x) myfun(x));
        pass(n, 5) = isnan(f);
    catch ME
        pass(n, 5) = strcmpi(ME.message, ...
            'Function returned NaN when evaluated.');
    end
    
    % Test a non-adaptive construction
    % [TODO]:  This test fails for funcheb1.
    p.funcheb.n = 11;
    try
        f = testclass.make(@(x) myfun(x), p);
        pass(n, 6) = isnan(f);
    catch ME
        pass(n, 6) = strcmpi(ME.message, ...
            'Function returned NaN when evaluated.');
    end

end

end


function y = myfun(x)
    y = 1./x;
    y(x == 0) = NaN;
end

% Test file for chebtech/isfinite.m

function pass = test_isfinite(pref)

if ( nargin < 1 )
    pref = chebtech.techPref();
end
p = pref;

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    % Test a scalar-valued function:
    y = ones(11, 1); 
    y(4) = inf; 
    f = testclass.make({[], y}, [], []);
    pass(n, 1) = ~isfinite(f);
    
    % Test an array-valued function:
    y = ones(11,1); 
    y(4) = inf; 
    f = testclass.make({[], y}, [], p);
    pass(n, 2) = ~isfinite(f);
    
    % Test a finite scalar-valued function:
    p = pref;
    f = testclass.make(@(x) x, [], p);
    pass(n, 3) = isfinite(f);
    
    % Test a finite array-valued function:
    f = testclass.make(@(x) [x, x], [], p);
    pass(n, 4) = isfinite(f);
end

end

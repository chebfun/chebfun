% Test file for chebtech/isreal.m

function pass = test_isreal(pref)

if ( nargin < 1 )
    pref = chebtech.techPref();
end

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    % Test a scalar-valued function:
    f = testclass.make(@(x) sin(x) + 1i*cos(x), [], pref);
    pass(n, 1) = ~isreal(f);
    
    f = testclass.make(@(x) 1i*cos(x), [], pref);
    pass(n, 2) = ~isreal(f);
    
    f = testclass.make(@(x) sin(x), [], pref);
    pass(n, 3) = isreal(f);
    
    % Test an array-valued function:
    f = testclass.make(@(x) [sin(x) + 1i*cos(x), exp(x)], [], pref);
    pass(n, 4) = ~isreal(f);
    
    f = testclass.make(@(x) [1i*cos(x), exp(x)], [], pref);
    pass(n, 5) = ~isreal(f);
    
    f = testclass.make(@(x) [sin(x), exp(x)], [], pref);
    pass(n, 6) = isreal(f);
end

end

% Test file for funcheb/isreal.m

function pass = test_isreal(pref)

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
    f = testclass.make(@(x) sin(x) + 1i*cos(x), [], [], p);
    pass(n, 1) = ~isreal(f);
    
    f = testclass.make(@(x) 1i*cos(x), [], [], p);
    pass(n, 2) = ~isreal(f);
    
    f = testclass.make(@(x) sin(x), [], [], p);
    pass(n, 3) = isreal(f);
    
    % Test a multi-valued function:
    f = testclass.make(@(x) [sin(x) + 1i*cos(x), exp(x)], [], [], p);
    pass(n, 4) = ~isreal(f);
    
    f = testclass.make(@(x) [1i*cos(x), exp(x)], [], [], p);
    pass(n, 5) = ~isreal(f);
    
    f = testclass.make(@(x) [sin(x), exp(x)], [], [], p);
    pass(n, 6) = isreal(f);
end

end

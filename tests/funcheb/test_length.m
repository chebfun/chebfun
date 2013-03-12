function pass = test_length(pref)

if ( nargin < 1 )
    pref = funcheb.pref;
end

for ( n = 1:2 )
    if ( n == 1 )
        testclass = funcheb1();
    else 
        testclass = funcheb2();
    end

    f = testclass.make(@(x) sin(x), 0, pref);
    pass(n, 1) = length(f) == 14;
    
    f = testclass.make(@(x) [sin(x), cos(x), 1i*exp(x)], 0, pref);
    pass(n, 2) = length(f) == 15;
    
    p = pref;
    p.funcheb.n = 101;
    f = testclass.make(@(x) [sin(x), cos(x), 1i*exp(x)], 0, p);
    pass(n, 3) = length(f) == 101;
end

end

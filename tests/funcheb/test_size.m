% Test file for funcheb/size.m

function pass = test_size(pref)

if ( nargin < 1 )
    pref = funcheb.pref;
end

for ( n = 1:2 )
    if ( n == 1 )
        testclass = funcheb1();
    else 
        testclass = funcheb2();
    end

    f = testclass.make(@(x) sin(x), [], [], pref);
    pass(n, 1) = all(size(f) == [14, 1]);
    
    f = testclass.make(@(x) [sin(x), cos(x), 1i*exp(x)], [], [], pref);
    pass(n, 2) = all(size(f) == [15, 3]);
    
    p = pref;
    p.funcheb.n = 101;
    f = testclass.make(@(x) [sin(x), cos(x), 1i*exp(x)], [], [], p);
    pass(n, 3) = all(size(f) == [101, 3]);
end

end

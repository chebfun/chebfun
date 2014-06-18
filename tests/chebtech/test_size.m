% Test file for chebtech/size.m

function pass = test_size(pref)

if ( nargin < 1 )
    pref = chebtech.techPref();
end

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    f = testclass.make(@(x) sin(x), [], pref);
    pass(n, 1) = all(size(f) == size(f.coeffs));
    
    f = testclass.make(@(x) [sin(x), cos(x), 1i*exp(x)], [], pref);
    pass(n, 2) = all(size(f) == size(f.coeffs));
    
    p = pref;
    p.fixedLength = 101;
    f = testclass.make(@(x) [sin(x), cos(x), 1i*exp(x)], [], p);
    pass(n, 3) = all(size(f) == [101, 3]);
end

end

% Test file for chebtech/length.m

function pass = test_length(pref)

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
    pass(n, 1) = length(f) == size(f.coeffs, 1);
    
    f = testclass.make(@(x) [sin(x), cos(x), 1i*exp(x)], [], pref);
    pass(n, 2) = length(f) == size(f.coeffs, 1);
    
    p = pref;
    p.fixedLength = 101;
    f = testclass.make(@(x) [sin(x), cos(x), 1i*exp(x)], [], p);
    pass(n, 3) = length(f) == 101;
end

end

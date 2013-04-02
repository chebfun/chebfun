% Test file for chebtech/length.m

function pass = test_length(pref)

if ( nargin < 1 )
    pref = chebtech.pref;
end

for ( n = 1:2 )
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    f = testclass.make(@(x) sin(x), [], [], pref);
    pass(n, 1) = length(f) == 14;
    
    f = testclass.make(@(x) [sin(x), cos(x), 1i*exp(x)], [], [], pref);
    pass(n, 2) = length(f) == 15;
    
    p = pref;
    p.chebtech.n = 101;
    f = testclass.make(@(x) [sin(x), cos(x), 1i*exp(x)], [], [], p);
    pass(n, 3) = length(f) == 101;
end

end

% Test file for fun/size.m

function pass = test_size(pref)

if ( nargin < 1 )
    pref = fun.pref;
end

pass = zeros(1, 3); % Pre-allocate pass matrix
for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
    else 
        testclass = unbndfun();
    end
    
    dom = [-2 7];
    f = testclass.make(@(x) sin(x), dom, [], [], pref);
    pass(n, 1) = all(size(f) == size(f.onefun.values));
    
    f = testclass.make(@(x) [sin(x), cos(x), 1i*exp(x)], dom, [], [], pref);
    pass(n, 2) = all(size(f) == size(f.onefun.values));
    
    pref = chebtech.pref(pref);
    pref.chebtech.n = 101;
    f = testclass.make(@(x) [sin(x), cos(x), 1i*exp(x)], dom, [], [], pref);
    pass(n, 3) = all(size(f) == [101, 3]);
end

end

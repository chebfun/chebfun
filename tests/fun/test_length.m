% Test file for fun/length.m

function pass = test_length(pref)

if ( nargin < 1 )
    pref = fun.pref;
end

pass = zeros(1, 3); % Pre-allocate pass matrix
for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
        dom = [-2 7];
    else 
        testclass = unbndfun();
    end
    
    f = testclass.make(@(x) sin(x), dom, [], [], pref);
    pass(n, 1) = length(f) == size(f.onefun.values, 1);
    
    f = testclass.make(@(x) [sin(x), cos(x), 1i*exp(x)], dom, [], [], pref);
    pass(n, 2) = length(f) == size(f.onefun.values, 1);
    
    pref = chebtech.pref(pref);
    pref.chebtech.n = 101;
    f = testclass.make(@(x) [sin(x), cos(x), 1i*exp(x)], dom, [], [], pref);
    pass(n, 3) = length(f) == 101;
end

end
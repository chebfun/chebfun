% Test file for fun/isreal.m

function pass = test_isreal(pref)

if ( nargin < 1 )
    pref = fun.pref;
end
p = pref;

pass = zeros(1, 6); % Pre-allocate pass matrix
for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
        dom = [-2 7];
    else
        testclass = unbndfun();
    end
    
    % Test a scalar-valued function:
    f = testclass.make(@(x) sin(x) + 1i*cos(x), dom, [], [], p);
    pass(n, 1) = ~isreal(f);
    
    f = testclass.make(@(x) 1i*cos(x), dom, [], [], p);
    pass(n, 2) = ~isreal(f);
    
    f = testclass.make(@(x) sin(x), dom, [], [], p);
    pass(n, 3) = isreal(f);
    
    % Test an array-valued function:
    f = testclass.make(@(x) [sin(x) + 1i*cos(x), exp(x)], dom, [], [], p);
    pass(n, 4) = ~isreal(f);
    
    f = testclass.make(@(x) [1i*cos(x), exp(x)], dom, [], [], p);
    pass(n, 5) = ~isreal(f);
    
    f = testclass.make(@(x) [sin(x), exp(x)], dom, [], [], p);
    pass(n, 6) = isreal(f);
end

end

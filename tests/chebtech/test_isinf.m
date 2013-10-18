% Test file for chebtech/isinf.m

function pass = test_isinf(pref)

if ( nargin < 1 )
    pref = chebtech.pref();
end
p = pref;

pass = zeros(2, 4); % Pre-allocate pass matrix
for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    % Test a scalar-valued function:
    p.numSamples = 11;
    y = testclass.chebpts(p.numSamples); % Force singularity to fall on grid.
    f = testclass.make(@(x) 1./(x - y(4)), [], [], p);
    pass(n, 1) = isinf(f);
    
    % Test an array-valued function:
    p.numSamples = 11;
    y = testclass.chebpts(p.numSamples); % Force singularity to fall on grid.
    f = testclass.make(@(x) [1./(x - y(4)), x], [], [], p);
    pass(n, 2) = isinf(f);
    
    % Test a finite scalar-valued function:
    p = pref;
    f = testclass.make(@(x) x, [], [], p);
    pass(n, 3) = ~isinf(f);
    
    % Test a finite array-valued function:
    f = testclass.make(@(x) [x, x], [], [], p);
    pass(n, 4) = ~isinf(f);
end

end

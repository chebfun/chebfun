% Test file for chebtech/isfinite.m

function pass = test_isfinite(pref)

if ( nargin < 1 )
    pref = chebtech.pref;
end
p = pref;

pass = zeros(2,4); % Pre-allocate pass matrix
for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    % Test a scalar-valued function:
    p.chebtech.n = 11;
    y = testclass.chebpts(p.chebtech.n); % Force singularity to fall on grid.
    f = testclass.make(@(x) 1./(x - y(4)), [], [], p);
    pass(n, 1) = ~isfinite(f);
    
    % Test a vector-valued function:
    p.chebtech.n = 11;
    y = testclass.chebpts(p.chebtech.n); % Force singularity to fall on grid.
    f = testclass.make(@(x) [1./(x - y(4)), x], [], [], p);
    pass(n, 2) = ~isfinite(f);
    
    % Test a finite scalar-valued function:
    p = pref;
    f = testclass.make(@(x) x, [], [], p);
    pass(n, 3) = isfinite(f);
    
    % Test a finite vector-valued function:
    f = testclass.make(@(x) [x, x], [], [], p);
    pass(n, 4) = isfinite(f);
end

end

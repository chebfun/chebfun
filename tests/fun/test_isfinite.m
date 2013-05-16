% Test file for fun/isfinite.m

function pass = test_isfinite(pref)

if ( nargin < 1 )
    pref = fun.pref;
end

pass = zeros(1, 4); % Pre-allocate pass matrix
for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();    
        dom = [-2 7];
    else 
        testclass = unbndfun();
    end

    % Test a scalar-valued function:
    pref = chebtech.pref(pref);
    pref.chebtech.n = 11;

    y = chebtech2.chebpts(pref.chebtech.n); % Force singularity to fall on grid.
    map = testclass.createMap(dom);
    y = map.for(y);
    f = testclass.make(@(x) 1./(x - y(4)), dom, [], [], pref);
    pass(n, 1) = ~isfinite(f);
    
    % Test a vector-valued function:
    f = testclass.make(@(x) [1./(x - y(4)), x], dom, [], [], pref);
    pass(n, 2) = ~isfinite(f);
    
    % Test a finite scalar-valued function:
    p = pref;
    f = testclass.make(@(x) x, dom, [], [], pref);
    pass(n, 3) = isfinite(f);
    
    % Test a finite vector-valued function:
    f = testclass.make(@(x) [x, x], dom, [], [], pref);
    pass(n, 4) = isfinite(f);
end

end
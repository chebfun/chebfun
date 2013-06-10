% Test file for fun/isnan.m

function pass = test_isnan(pref)

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
    f = testclass.make(@(x) x, dom, [], [], p);
    pass(n, 1) = ~isnan(f);
    
    % Test a vector-valued function:
    f = testclass.make(@(x) [x, x.^2], dom, [], [], p);
    pass(n, 2) = ~isnan(f);
    
    % Test a NaN scalar-valued function:
    try
        f = testclass.make(@(x) x + NaN, dom, [], [], p);
        pass(n, 3) = isnan(f);
    catch ME
        pass(n, 3) = strcmpi(ME.message, 'Too many NaNs/Infs to handle.');
    end
    
    % Test a NaN vector-valued function:
    try
        f = testclass.make(@(x) [x, x + NaN], dom, [], [], p);
        pass(n, 4) = isnan(f);
    catch ME
        pass(n, 4) = strcmpi(ME.message, 'Too many NaNs/Infs to handle.');
    end

    % Test a NaN vector-valued function:
    try
        f = testclass.make(@(x) myfun(x), dom, [], [], p);
        pass(n, 5) = isnan(f);
    catch ME
        pass(n, 5) = strcmpi(ME.message, ...
            'Function returned NaN when evaluated.');
    end

    % Test a non-adaptive construction
    p = chebtech.pref(p);
    p.chebtech.n = 11;
    try
        f = testclass.make(@(x) myfun(x), dom, [], [], p);
        pass(n, 6) = isnan(f);
    catch ME
        pass(n, 6) = strcmpi(ME.message, ...
            'Function returned NaN when evaluated.');
    end

end

end


function y = myfun(x)
    y = 1./x;

    % Put NaN at point closest to zero on the sampling grid.  Can't just
    % put it at zero directly, since that won't work for chebtech1.
    [ignored, n] = min(abs(x));
    y(n) = NaN;
end

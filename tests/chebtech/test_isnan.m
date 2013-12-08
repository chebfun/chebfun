% Test file for chebtech/isnan.m

function pass = test_isnan(pref)

if ( nargin < 1 )
    pref = chebtech.techPref();
end
p = pref;

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    % Test a scalar-valued function:
    f = testclass.make(@(x) x, [], [], p);
    pass(n, 1) = ~isnan(f);

    % Test an array-valued function:
    f = testclass.make(@(x) [x, x.^2], [], [], p);
    pass(n, 2) = ~isnan(f);

    % Artificially construct and test a NaN-valued function:
    f = testclass.make(NaN);
    pass(n, 3) = isnan(f);

    % Test a NaN scalar-valued function:
    try
        f = testclass.make(@(x) x + NaN);
        pass(n, 4) = isnan(f);
    catch ME
        pass(n, 4) = strcmpi(ME.message, 'Too many NaNs/Infs to handle.');
    end

    % Test a NaN array-valued function:
    try
        f = testclass.make(@(x) [x, x + NaN]);
        pass(n, 5) = isnan(f);
    catch ME
        pass(n, 5) = strcmpi(ME.message, 'Too many NaNs/Infs to handle.');
    end

    % Test a NaN array-valued function:
    try
        f = testclass.make(@(x) myfun(x));
        pass(n, 6) = isnan(f);
    catch ME
        pass(n, 6) = strcmpi(ME.message, ...
            'Function returned NaN when evaluated.');
    end

    % Test a non-adaptive construction
    p.numPoints = 11;
    try
        f = testclass.make(@(x) myfun(x), [], [], p);
        pass(n, 7) = isnan(f);
    catch ME
        pass(n, 7) = strcmpi(ME.message, ...
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

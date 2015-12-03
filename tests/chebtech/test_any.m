% Test file for @chebtech/any.m.

function pass = test_any(pref)

if ( nargin < 1 )
    pref = chebtech.techPref();
end

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else
        testclass = chebtech2();
    end

    % Check behavior for any() down columns.
    pass(n, 1) = ~any(testclass);

    f = testclass.make(@(x) [sin(x) 0*x cos(x)]);
    pass(n, 2) = isequal(any(f), [1 0 1]);

    % Check behavior for any() across rows.
    g = any(f, 2);
    pass(n, 3) = isequal(g.coeffs, 1);

    f = testclass.make(@(x) [0*x 0*x]);
    g = any(f, 2);
    pass(n, 4) = isequal(g.coeffs, 0);
end

end

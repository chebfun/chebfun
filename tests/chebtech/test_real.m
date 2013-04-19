% Test file for chebtech/real.m

function pass = test_real(pref)

if ( nargin < 1 )
    pref = chebtech.pref;
end
tol = 10*pref.chebtech.eps;

pass = zeros(2,4); % Pre-allocate pass matrix
for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    % Test a scalar-valued function:
    f = testclass.make(@(x) cos(x) + 1i*sin(x), [], [], pref);
    g = testclass.make(@(x) cos(x), [], [], pref);
    h = real(f);
    pass(n, 1) = norm(h.values - g.values, inf) < tol;
    
    % Test a multi-valued function:
    f = testclass.make(@(x) [cos(x) + 1i*sin(x), -exp(1i*x)], [], [], pref);
    g = testclass.make(@(x) [cos(x), -real(exp(1i*x))], [], [], pref);
    h = real(f);
    pass(n, 2) = norm(h.values - g.values, inf) < tol;
    
    % Test a real function:
    f = testclass.make(@(x) 1i*cos(x), [], [], pref);
    g = real(f);
    pass(n, 3) = numel(g.values) == 1 && g.values == 0;
    
    
    % Test a multivalued real function:
    f = testclass.make(@(x) 1i*[cos(x), sin(x), exp(x)], [], [], pref);
    g = real(f);
    pass(n, 4) = all(size(g.values) == [1, 3]) && all(g.values == 0);
end

end

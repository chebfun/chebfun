% Test file for fun/real.m

function pass = test_real(pref)

if ( nargin < 1 )
    pref = fun.pref;
end
tol = 10*pref.fun.eps;

pass = zeros(1, 4); % Pre-allocate pass matrix
for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
        dom = [-2 7];
    else 
        testclass = unbndfun();
    end

    % Test a scalar-valued function:
    f = testclass.make(@(x) exp(1i*x) + 1i*sin(x), dom, [], [], pref);
    g = testclass.make(@(x) cos(x), dom, [], [], pref);
    h = real(f);
    pass(n, 1) = norm(h.onefun.values - g.onefun.values, inf) < tol;
    
    % Test an array-valued function:
    f = testclass.make(@(x) [exp(1i*x) + 1i*sin(x), -exp(1i*x)], dom, [], [], pref);
    g = testclass.make(@(x) [cos(x), -real(exp(1i*x))], dom, [], [], pref);
    h = real(f);
    pass(n, 2) = norm(h.onefun.values - g.onefun.values, inf) < tol;
    
    % Test a real function:
    f = testclass.make(@(x) 1i*cos(x), dom, [], [], pref);
    g = real(f);
    pass(n, 3) = numel(g.onefun.values) == 1 && g.onefun.values == 0;

    % Test an array-valued real function:
    f = testclass.make(@(x) 1i*[cos(x), sin(x), exp(x)], dom, [], [], pref);
    g = real(f);
    pass(n, 4) = all(size(g.onefun.values) == [1, 3]) && all(g.onefun.values == 0);
end

end

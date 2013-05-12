% Test file for fun/conj.

function pass = test_conj(pref)

if ( nargin < 1 )
    pref = fun.pref;
end

pass = zeros(1,2);

for n = 1:1 % [TODO]: unbndfun
    if (n == 1)
        testclass = bndfun();
    else 
        testclass = unbndfun();
    end

    tol = 10*pref.fun.eps;
    dom = [-2 7];
    
    % Test a scalar-valued function:
    f = testclass.make(@(x) cos(x) + 1i*sin(x), dom, [], [], pref);
    g = testclass.make(@(x) cos(x) - 1i*sin(x), dom, [], [], pref);
    h = conj(f);
    pass(n, 1) = norm(h.onefun.values - g.onefun.values, inf) < tol;
    
    % Test an array-valued function:
    f = testclass.make(@(x) [cos(x) + 1i*sin(x), -exp(1i*x)], dom, [], [], pref);
    g = testclass.make(@(x) cos(x) - 1i*sin(x), dom, [], [], pref);
    h = conj(f);
    pass(n, 2) = norm(h.onefun.values - [g.onefun.values, -(g.onefun.values)], inf) < tol;
end

end
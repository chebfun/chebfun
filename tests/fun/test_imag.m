% Test file for fun/imag.m

function pass = test_imag(pref)

if ( nargin < 1 )
    pref = fun.pref;
end

tol = 500*pref.fun.eps;

pass = zeros(1, 4); % Pre-allocate pass matrix
for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();    
        dom = [-2 7];
    else 
        testclass = unbndfun();
    end

    % Test a scalar-valued function:
    f = testclass.make(@(x) exp(x) + 1i*sin(x), dom, [], [], pref);
    g = testclass.make(@(x) sin(x), dom, [], [], pref);
    h = imag(f);
    h.onefun = prolong(h.onefun, length(g));
    pass(n, 1) = norm(h.onefun.values - g.onefun.values, inf) < tol;
    
    % Test an array-valued function:
    f = testclass.make(@(x) [exp(x) + 1i*sin(x), -exp(1i*x)], dom, [], [], pref);
    g = testclass.make(@(x) [sin(x), -imag(exp(1i*x))], dom, [], [], pref);
    h = imag(f);
    h.onefun = prolong(h.onefun, length(g));
    pass(n, 2) = norm(h.onefun.values - g.onefun.values, inf) < tol;
    
    % Test a real function:
    f = testclass.make(@(x) cos(x), dom, [], [], pref);
    g = imag(f);
    pass(n, 3) = numel(g.onefun.values) == 1 && g.onefun.values == 0;
    
    % Test an array-valued real function:
    f = testclass.make(@(x) [cos(x), sin(x), exp(x)], dom, [], [], pref);
    g = imag(f);
    pass(n, 4) = all(size(g.onefun.values) == [1, 3]) && all(g.onefun.values == 0);
end

end
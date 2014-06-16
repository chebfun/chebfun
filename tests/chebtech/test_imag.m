% Test file for chebtech/imag.m

function pass = test_imag(pref)

if ( nargin < 1 )
    pref = chebtech.techPref();
end

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    % Test a scalar-valued function:
    f = testclass.make(@(x) exp(x) + 1i*sin(x), [], pref);
    g = testclass.make(@(x) sin(x), [], pref);
    h = imag(f);
    h = prolong(h, length(g));
    pass(n, 1) = norm(h.coeffs - g.coeffs, inf) < 10*h.vscale.*h.epslevel;
    
    % Test an array-valued function:
    f = testclass.make(@(x) [exp(x) + 1i*sin(x), -exp(1i*x)], [], pref);
    g = testclass.make(@(x) [sin(x), -imag(exp(1i*x))], [], pref);
    h = imag(f);
    h = prolong(h, length(g));
    pass(n, 2) = norm(h.coeffs - g.coeffs, inf) < 10*max(h.vscale.*h.epslevel);
    
    % Test a real function:
    f = testclass.make(@(x) cos(x), [], pref);
    g = imag(f);
    pass(n, 3) = numel(g.coeffs) == 1;
    
    % Test an array-valued real function:
    f = testclass.make(@(x) [cos(x), sin(x), exp(x)], [], pref);
    g = imag(f);
    pass(n, 4) = all(size(g.coeffs) == [1, 3]) && all(g.coeffs == 0);
end

end

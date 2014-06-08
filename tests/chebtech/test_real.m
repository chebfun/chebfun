% Test file for chebtech/real.m

function pass = test_real(pref)

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
    f = testclass.make(@(x) exp(1i*x) + 1i*sin(x), [], pref);
    g = testclass.make(@(x) cos(x), [], pref);
    h = real(f);
    pass(n, 1) = norm(h.coeffs - g.coeffs, inf) < 10*h.vscale.*h.epslevel;
    
    % Test an array-valued function:
    f = testclass.make(@(x) [exp(1i*x) + 1i*sin(x), -exp(1i*x)], [], pref);
    g = testclass.make(@(x) [cos(x), -real(exp(1i*x))], [], pref);
    h = real(f);
    pass(n, 2) = norm(h.coeffs - g.coeffs, inf) < 10*max(h.vscale.*h.epslevel);
    
    % Test a real function:
    f = testclass.make(@(x) 1i*cos(x), [], pref);
    g = real(f);
    pass(n, 3) = numel(g.coeffs) == 1 && g.coeffs == 0;

    % Test an array-valued real function:
    f = testclass.make(@(x) 1i*[cos(x), sin(x), exp(x)], [], pref);
    g = real(f);
    pass(n, 4) = all(size(g.coeffs) == [1, 3]) && all(g.coeffs == 0);
end

end

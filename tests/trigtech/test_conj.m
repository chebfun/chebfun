% Test file for trigtech/conj.m

function pass = test_conj(pref)

if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

% Test a scalar-valued function:
f = testclass.make(@(x) cos(pi*x) + 1i*sin(pi*x), [], pref);
g = testclass.make(@(x) cos(pi*x) - 1i*sin(pi*x), [], pref);
h = conj(f);
pass(1) = norm(h.coeffs - g.coeffs, inf) < 10*vscale(h)*eps;

% Test an array-valued function:
f = testclass.make(@(x) [cos(pi*x) + 1i*sin(pi*x), -exp(1i*pi*x)], [], pref);
g = testclass.make(@(x) cos(pi*x) - 1i*sin(pi*x), [], pref);
h = conj(f);
pass(2) = norm(h.coeffs - [g.coeffs, -(g.coeffs)], inf) < ...
    10*max(vscale(h)*eps);

% Test an array-valued function with a real column and complex column
f = testclass.make(@(x) [exp(cos(pi*x)) exp(1i*pi*x)], [], pref);
g = testclass.make(@(x) [exp(cos(pi*x)) cos(pi*x) - 1i*sin(pi*x)], [], pref);
h = conj(f);
pass(3) = norm(h.coeffs - g.coeffs, inf) < ...
    10*max(vscale(h)*eps);

end

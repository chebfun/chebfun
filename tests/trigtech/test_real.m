% Test file for trigtech/real.m

function pass = test_real(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

% Test a scalar-valued function.
f = testclass.make(@(x) exp(20i*pi*x) + 1i*sin(100*pi*x), [], pref);
g = testclass.make(@(x) cos(20*pi*x), [], pref);
h = real(f);
g = prolong(g, length(h));
pass(1) = norm(h.coeffs - g.coeffs, inf) < 10*h.vscale.*h.epslevel;

% Test an array-valued function.
f = testclass.make(@(x) [exp(20i*pi*x) + 1i*sin(100*pi*x), -exp(10i*pi*x)], [], pref);
g = testclass.make(@(x) [cos(20*pi*x), -real(exp(10i*pi*x))], [], pref);
h = real(f);
pass(2) = norm(h.coeffs - g.coeffs, inf) < 10*max(h.vscale.*h.epslevel);

% Test a real function.
f = 1i*testclass.make(@(x) cos(30*pi*x), [], pref);
g = real(f);
pass(3) = numel(g.coeffs) == 1 && g.coeffs == 0;

% Test an array-valued real function.
f = 1i*testclass.make(@(x) [cos(99*pi*x), sin(99*pi*x), exp(cos(pi*x))], [], pref);
g = real(f);
pass(4) = all(size(g.coeffs) == [1, 3]) && all(g.coeffs == 0);
    
end

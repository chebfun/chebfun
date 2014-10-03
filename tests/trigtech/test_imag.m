% Test file for trigtech/imag.m

function pass = test_imag(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

%%
% Test a scalar-valued function.
f = testclass.make(@(x) cos(pi*x) + 1i*sin(pi*x), [], pref);
g = testclass.make(@(x) sin(pi*x), [], pref);
h = imag(f);
g = prolong(g, length(h));
pass(1) = norm(h.coeffs - g.coeffs, inf) < 10*h.vscale.*h.epslevel;

%%
% Test an array-valued function.
f = testclass.make(@(x) [cos(sin(pi*x)) + 1i*sin(cos(pi*x)), -exp(1i*pi*x)], [], pref);
g = testclass.make(@(x) [sin(cos(pi*x)), -imag(exp(1i*pi*x))], [], pref);
h = imag(f);
g = prolong(g, length(h));
pass(2) = norm(h.coeffs - g.coeffs, inf) < 100*max(h.vscale.*h.epslevel);

%%
% Test a real function.
f = testclass.make(@(x) cos(pi*x), [], pref);
g = imag(f);
pass(3) = numel(g.coeffs) == 1;

%%
% Test an array-valued real function.
f = testclass.make(@(x) [cos(pi*x), sin(pi*x)], [], pref);
g = imag(f);
pass(4) = all(size(g.coeffs) == [1, 2]) && all(g.coeffs == 0);

end

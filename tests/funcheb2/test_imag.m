function pass = test_imag(pref)

if ( nargin < 1 )
    pref = funcheb.pref;
end
tol = 10*pref.funcheb.eps;

% Test a scalar-valued function:
f = funcheb2(@(x) cos(x) + 1i*sin(x), pref);
g = funcheb2(@(x) sin(x), 0, pref);
h = imag(f);
h = prolong(h, length(g));
pass(1) = norm(h.values - g.values, inf) < tol;

% Test a multi-valued function:
f = funcheb2(@(x) [cos(x) + 1i*sin(x), -exp(1i*x)], pref);
g = funcheb2(@(x) [sin(x), -imag(exp(1i*x))], 0, pref);
h = imag(f);
h = prolong(h, length(g));
pass(2) = norm(h.values - g.values, inf) < tol;

% Test a real function:
f = funcheb2(@(x) cos(x), pref);
g = imag(f);
pass(3) = numel(g.values) == 1 && g.values == 0;

% Test a multivalued real function:
f = funcheb2(@(x) [cos(x), sin(x), exp(x)], pref);
g = imag(f);
pass(4) = all(size(g.values) == [1, 3]) && all(g.values == 0);

end

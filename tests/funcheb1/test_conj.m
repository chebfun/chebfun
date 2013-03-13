function pass = test_conj(pref)

if ( nargin < 1 )
    pref = funcheb.pref;
end
tol = 10*pref.funcheb.eps;

% Test a scalar-valued function:
f = funcheb1(@(x) cos(x) + 1i*sin(x), 0, pref);
g = funcheb1(@(x) cos(x) - 1i*sin(x), 0, pref);
h = conj(f);
pass(1) = norm(h.values - g.values, inf) < tol;

% Test a multi-valued function:
f = funcheb1(@(x) [cos(x) + 1i*sin(x), -exp(1i*x)], 0, pref);
g = funcheb1(@(x) cos(x) - 1i*sin(x), 0, pref);
h = conj(f);
pass(2) = norm(h.values - [g.values, -(g.values)], inf) < tol;

end

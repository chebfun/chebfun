function pass = test_real(pref)

if ( nargin < 1 )
    pref = funcheb.pref;
end
tol = 10*pref.funcheb.eps;

% Test a scalar-valued function:
f = funcheb1(@(x) cos(x) + 1i*sin(x), 0, pref);
g = funcheb1(@(x) cos(x), 0, pref);
h = real(f);
pass(1) = norm(h.values - g.values, inf) < tol;

% Test a multi-valued function:
f = funcheb1(@(x) [cos(x) + 1i*sin(x), -exp(1i*x)], 0, pref);
g = funcheb1(@(x) [cos(x), -real(exp(1i*x))], 0, pref);
h = real(f);
pass(2) = norm(h.values - g.values, inf) < tol;

% Test a real function:
f = funcheb1(@(x) 1i*cos(x), pref);
g = real(f);
pass(3) = numel(g.values) == 1 && g.values == 0;


% Test a multivalued real function:
f = funcheb1(@(x) 1i*[cos(x), sin(x), exp(x)], pref);
g = real(f);
pass(4) = all(size(g.values) == [1, 3]) && all(g.values == 0);

end

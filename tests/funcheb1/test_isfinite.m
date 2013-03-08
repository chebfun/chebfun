function pass = test_isfinite(pref)

if ( nargin < 1 )
    pref = funcheb1.pref;
end
p = pref;

% Test a scalar-valued function:
p.funcheb1.n = 11;
y = funcheb1.chebpts(p.funcheb1.n);
f = funcheb1(@(x) 1./(x-y(4)), 0, p);
pass(1) = ~isfinite(f);

% Test a vector-valued function:
p.funcheb1.n = 11;
f = funcheb1(@(x) [1./(x-y(4)), x], 0, p);
pass(2) = ~isfinite(f);

% Test a finite scalar-valued function:
p = pref;
f = funcheb1(@(x) x, 0, p);
pass(3) = isfinite(f);

% Test a finite vector-valued function:
f = funcheb1(@(x) [x, x], 0, p);
pass(4) = isfinite(f);

end

function pass = test_isreal(pref)

if ( nargin < 1 )
    pref = funcheb1.pref;
end
p = pref;

% Test a scalar-valued function:
f = funcheb1(@(x) sin(x) + 1i*cos(x), 0, p);
pass(1) = ~isreal(f);

f = funcheb1(@(x) 1i*cos(x), 0, p);
pass(2) = ~isreal(f);

f = funcheb1(@(x) sin(x), 0, p);
pass(3) = isreal(f);

% Test a multi-valued function:
f = funcheb1(@(x) [sin(x) + 1i*cos(x), exp(x)], 0, p);
pass(4) = ~isreal(f);

f = funcheb1(@(x) [1i*cos(x), exp(x)], 0, p);
pass(5) = ~isreal(f);

f = funcheb1(@(x) [sin(x), exp(x)], 0, p);
pass(6) = isreal(f);


end

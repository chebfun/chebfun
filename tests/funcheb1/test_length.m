function pass = test_length(pref)

if ( nargin < 1 )
    pref = funcheb.pref;
end

f = funcheb1(@(x) sin(x), 0, pref);
pass(1) = length(f) == 14;

f = funcheb1(@(x) [sin(x), cos(x), 1i*exp(x)], 0, pref);
pass(2) = length(f) == 15;

p = pref;
p.funcheb.n = 101;
f = funcheb1(@(x) [sin(x), cos(x), 1i*exp(x)], 0, p);
pass(3) = length(f) == 101;


end
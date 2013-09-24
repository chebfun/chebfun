% Test file for @chebfun/all.m.

function pass = test_all(pref)

if (nargin < 1)
    pref = chebfun.pref();
end

% Test scalar valued chebfuns.
f = chebfun(@(x) sin(x), [-1 -0.5 0 0.5 1], pref);
pass(1) = ~all(f);

f = chebfun(@(x) sin(x - 0.1), [-1 -0.5 0 0.5 1], pref);
pass(2) = ~all(f);

f = chebfun(@(x) exp(2*pi*1i*x), [-1 -0.5 0 0.5 1], pref);
pass(3) = all(f);

% Test array-valued chebfun.
f = chebfun(@(x) [sin(x) sin(x - 0.1) exp(2*pi*1i*x)], [-1 -0.5 0 0.5 1], pref);
pass(4) = isequal(all(f), logical([0 0 1]));

end

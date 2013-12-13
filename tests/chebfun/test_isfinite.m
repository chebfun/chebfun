% Test file for @chebfun/isfinite.m.

function pass = test_isfinite(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

% Test finite scalar function.
f = chebfun(@(x) sin(x));
pass(1) = isfinite(f);

% Test finite array-valued function.
f = chebfun(@(x) [sin(x) cos(x)]);
pass(2) = isfinite(f);

% Test on function with infinite breakpoint value.
% [TODO]:  Test with a less artificial example.
val = f.pointValues(1, 1);
f.pointValues(1, 1) = Inf;
pass(3) = ~isfinite(f);

% [TODO]:  Add a test with a singular function once we have singfun.

end

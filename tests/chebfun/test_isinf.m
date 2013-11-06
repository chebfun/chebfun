% Test file for @chebfun/isinf.m.

function pass = test_isinf(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

% Test finite scalar function.
f = chebfun(@(x) sin(x));
pass(1) = ~isinf(f);

% Test finite array-valued function.
f = chebfun(@(x) [sin(x) cos(x)]);
pass(2) = ~isinf(f);

% Test on function with infinite breakpoint value.
% [TODO]:  Test with a less artificial example.
val = f.impulses(1, 1);
f.impulses(1, 1) = Inf;
pass(3) = isinf(f);

% Test on function with nontrivial higher-order impulses.
% [TODO]:  Test with a less artificial example.
f.impulses(1, 1) = val;
f.impulses(:,:,2) = ones(size(f.impulses, 1), size(f.impulses, 2));
pass(4) = isinf(f);

% [TODO]:  Add a test with a singular function once we have singfun.

end

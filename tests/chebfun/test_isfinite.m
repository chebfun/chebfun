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
val = f.impulses(1, 1);
f.impulses(1, 1) = Inf;
pass(3) = ~isfinite(f);

% Test on function with nontrivial higher-order impulses.
% [TODO]:  Test with a less artificial example.
f.impulses(1, 1) = val;
f.impulses(:,:,2) = ones(size(f.impulses, 1), size(f.impulses, 2));
pass(4) = ~isfinite(f);

% Integration of SINGFUN:
dom = [-2 7];
pow = -1.64;
f = chebfun(@(x) sin(100*x).*(x-dom(1)).^pow, dom, 'exps', [pow 0], 'splitting', 'on');
pass(5) = ~isfinite(f);

end

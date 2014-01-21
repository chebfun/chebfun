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
f.pointValues(1, 1) = Inf;
pass(3) = isinf(f);


%% Test on singular function:
dom = [-2 7];
pow = -1.64;
f = chebfun(@(x) sin(100*x).*(x-dom(1)).^pow, dom, 'exps', [pow 0], ...
    'splitting', 'on');
pass(5) = isinf(f);

end
% Test file for @chebfun/isfinite.m.

function pass = test_isfinite(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Test finite scalar function.
f = chebfun(@(x) sin(x));
pass(1) = isfinite(f);

% Test finite array-valued function.
f = chebfun(@(x) [sin(x) cos(x)]);
pass(2) = isfinite(f);

% Test on function with infinite breakpoint value.
% [TODO]:  Test with a less artificial example.
f.pointValues(1, 1) = Inf;
pass(3) = ~isfinite(f);

%% Test on singular function:
dom = [-2 7];
pow = -1.64;
f = chebfun(@(x) sin(100*x).*(x-dom(1)).^pow, dom, 'exps', [pow 0], ...
    'splitting', 'on');
pass(4) = ~isfinite(f);

%% Test for functions defined on unbounded domain:

% Function defined on [0 Inf]:

% Specify the domain:
dom = [0 Inf];

op = @(x) 0.75+sin(10*x)./exp(x);
f = chebfun(op, dom, 'splitting', 'on');
pass(5) = isfinite(f);

% Function defined on [0 Inf]:

% Set the domain:
dom = [-Inf -3*pi];

% Blow-up function:
op = @(x) x.*(5+exp(x.^3))./(dom(2)-x);
f = chebfun(op, dom, 'exps', [0 -1]); 
pass(6) = ~isfinite(f);

end

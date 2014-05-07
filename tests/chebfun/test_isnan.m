% Test file for @chebfun/isnan.m.

function pass = test_isnan(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Check empty case.
pass(1) = ~isnan(chebfun());

% Check clearly non-NaN cases.
f = chebfun(@(x) sin(x), [-1 -0.5 0 0.5 1], pref);
pass(2) = ~isnan(f);
g = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);
pass(3) = ~isnan(g);

% Check NaN piontValues values.
f.pointValues(2,1) = NaN;
pass(4) = isnan(f);

% Check a case with a NaN fun.  This is an artificial construction, but it's
% the only way to do this at the moment.
nanfun = chebtech2(NaN);
f.funs{2}.onefun = nanfun;
pass(5) = isnan(f);

%% Test on singular function:
dom = [-2 7];
pow = -1.64;
f = chebfun(@(x) sin(100*x).*(x-dom(1)).^pow, dom, 'exps', [pow 0], ...
    'splitting', 'on');
pass(6) = ~isnan(f);

%% Test for function defined on unbounded domain:

% Function defined on [0 Inf]:

% Specify the domain:
dom = [0 Inf];

op = @(x) 0.75+sin(10*x)./exp(x);
f = chebfun(op, dom, 'splitting', 'on');
pass(7) = ~isnan(f);

% Function defined on [0 Inf]:

% Set the domain:
dom = [-Inf -3*pi];

% Blow-up function:
op = @(x) x.*(5+exp(x.^3))./(dom(2)-x);
f = chebfun(op, dom, 'exps', [0 -1]); 
pass(8) = ~isnan(f);

end

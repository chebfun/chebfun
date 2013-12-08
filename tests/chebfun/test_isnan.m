% Test file for @chebfun/isnan.m.

function pass = test_isnan(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

% Check empty case.
pass(1) = ~isnan(chebfun());

% Check clearly non-NaN cases.
f = chebfun(@(x) sin(x), [-1 -0.5 0 0.5 1], pref);
pass(2) = ~isnan(f);
g = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);
pass(3) = ~isnan(g);

% Check NaN impulse values.
f.impulses(2,1,1) = NaN;
pass(4) = isnan(f);
g.impulses(4,3,3) = NaN;
pass(5) = isnan(g);

% Check a case with a NaN fun.  This is an artificial construction, but it's
% the only way to do this at the moment.
nanfun = chebtech2(NaN);
f.funs{2}.onefun = nanfun;
pass(6) = isnan(f);

% Integration of SINGFUN:
dom = [-2 7];
pow = -1.64;
f = chebfun(@(x) sin(100*x).*(x-dom(1)).^pow, dom, 'exps', [pow 0], 'splitting', 'on');
pass(7) = ~isnan(f);

end
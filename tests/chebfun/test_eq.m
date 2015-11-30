% Test file for @chebfun/eq.m.

function pass = test_eq(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Check the empty case.
f = chebfun(@(x) sin(x), [-1 -0.5 0 0.5 1], pref);
g = chebfun();
pass(1) = isempty(f == g)  && isempty(g == f);

% Check a few simple examples.
g = chebfun(@(x) 0*x + sqrt(2)/2, pref);
h = f == g;

ind = find(h.pointValues == 1);
pass(2) = abs(h.domain(ind) - pi/4) < 10*eps;

f = chebfun(@(x) exp(x), pref);
g = chebfun(@(x) (exp(0.5) - exp(-0.5))*(x + 0.5) + exp(-0.5), pref);
h = f == g;

ind = find(h.pointValues == 1);
pass(3) = norm(h.domain(ind) - [-0.5 0.5], inf) < 10*eps;

h = f == f;
pass(4) = (numel(h.funs) == 1) && all(feval(h, x) == 1);

h = f == -f;
pass(5) = (numel(h.funs) == 1) && all(feval(h, x) == 0);

%% Check an example where eq() does not need to introduce any new breakpoints.
f = chebfun(@(x) exp(x), [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) (exp(0.5) - exp(-0.5))*(x + 0.5) + exp(-0.5), pref);
h = f == g;

ind = find(h.pointValues == 1);
pass(6) = norm(h.domain(ind) - [-0.5 0.5], inf) < 10*eps;

% Check error conditions.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) exp(2*pi*1i*x), [-1 1], pref);

try
    h = f == g
    pass(7) = false;
catch ME
    pass(7) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:eq:array');
end

try
    h = g == f
    pass(8) = false;
catch ME
    pass(8) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:eq:array');
end

% Test for singular function:

dom = [-2 7];

% Generate a few random points to use as test values:
x = diff(dom) * rand(100, 1) + dom(1);

f = chebfun(@(x) cos(10*x)./(x - dom(1)).^0.6, dom, 'exps', [-0.6 0]);
h = (f == f);
hVals = feval(h, x);
pass(8) = ~any(hVals - 1);

end

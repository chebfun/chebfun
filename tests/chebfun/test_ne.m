% Test file for @chebfun/ne.m.

function pass = test_ne(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Check the empty case.
f = chebfun(@(x) sin(x), [-1 -0.5 0 0.5 1], pref);
g = chebfun();
pass(1) = isempty(f ~= g)  && isempty(g ~= f);

% Check a few simple examples.
g = chebfun(@(x) 0*x + sqrt(2)/2, pref);
h = f ~= g;
ind = find(h.pointValues == 0);
pass(2) = abs(h.domain(ind) - pi/4) < 10*vscale(h)*eps;

f = chebfun(@(x) exp(x), pref);
g = chebfun(@(x) (exp(0.5) - exp(-0.5))*(x + 0.5) + exp(-0.5), pref);
h = f ~= g;
ind = find(h.pointValues == 0);
pass(3) = norm(h.domain(ind) - [-0.5 0.5], inf) < 10*vscale(h)*eps;

h = f ~= f;
pass(4) = (numel(h.funs) == 1) && all(feval(h, x) == 0);

h = f ~= -f;
pass(5) = (numel(h.funs) == 1) && all(feval(h, x) == 1);

%% Check an example where ne() does not need to introduce any new breakpoints.
f = chebfun(@(x) exp(x), [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) (exp(0.5) - exp(-0.5))*(x + 0.5) + exp(-0.5), pref);
h = f ~= g;
ind = find(h.pointValues == 0);
pass(6) = norm(h.domain(ind) - [-0.5 0.5], inf) < 10*vscale(h)*eps;

% Check error conditions.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) exp(2*pi*1i*x), [-1 1], pref);

try
    h = f ~= g
    pass(7) = false;
catch ME
    pass(7) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:ne:array');
end

try
    h = g ~= f
    pass(8) = false;
catch ME
    pass(8) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:ne:array');
end

%% Test for singular function:

f = chebfun(@(x) -sin(x)./(x+1), 'exps', [-1 0]);
g = ( f ~= f );
pass(9) = ( g.pointValues(1) == 0 && ~any(feval(g, x)) );

%% Test for function defined on unbounded domain:

% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];
domCheck = [-1e2 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% Blow-up function:
op = @(x) x.^2.*(1-exp(-x.^2));
f = chebfun(op, dom, 'exps', [2 2]);
g = ( f ~= f );
gVals = feval(g, x);
pass(10) = ~any(gVals);

end

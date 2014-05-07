% Test file for @chebfun/logical.m.

function pass = test_logical(pref)

if (nargin < 1)
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Check a few basic cases.
f = chebfun(@(x) sin(2*x), [-1 1], pref);
g = logical(f);
pass(1) = isequal(g.pointValues, [1 0 1].') && all(feval(g, x) == 1);

f = chebfun(@(x) exp(x), [-1 -0.5 0.5 1], pref);
g = logical(f);
pass(2) = all(g.pointValues == 1) && all(feval(g, x) == 1);

f = chebfun(@(x) 0*x, [-1 -0.5 0.5 1], pref);
g = logical(f);
pass(3) = all(g.pointValues == 0) && all(feval(g, x) == 0);

f = chebfun(@(x) 0*x + 1, [-1 -0.5 0.5 1], pref);
g = logical(f);
pass(4) = all(g.pointValues == 1) && all(feval(g, x) == 1);

% Check complex chebfun.
f = chebfun(@(x) exp(1i*x).*sin(x), [-1 -0.5 0.5 1], pref);
g = logical(f);
pass(5) = isequal(g.pointValues, [1 1 0 1 1].') && all(feval(g, x) == 1);

% Check array-valued chebfun.
f = chebfun(@(x) [sin(x) exp(x)], [-1 -0.5 0.5 1], pref);
g = logical(f);
pass(6) = isequal(g.pointValues, [1 1 0 1 1 ; 1 1 1 1 1].') && ...
    all(all(feval(g, x) == 1));

%% Test for singular function:

% define the domain:
dom = [-2 7];

op = @(x) sin(30*x)./((x-dom(1)).*(x-dom(2)));
f = chebfun(op, dom, 'exps', [-1 -1], 'splitting', 'on');
h = logical(f);

% check values:

% Generate a few random points to use as test values:
x = diff(dom) * rand(100, 1) + dom(1);

fval = feval(h, x);
err = fval - 1;
pass(7) = ~any( err );

r = roots(f);
pass(8) = ~any( h(r) );

%% Test for function defined on unbounded domain:

% Functions on [-inf b]:

% Set the domain:
dom = [-Inf 0 3*pi ];
domCheck = [-1e6 3*pi];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% Blow-up function:
op = @(x) 5+exp(x.^3);
opExact = @(x) logical(op(x));
f = chebfun({op 0}, dom);
g = logical(f);
gVals = feval(g, x);
gExact = opExact(x);
err = gVals - gExact;
pass(9) = ~any(err);

end

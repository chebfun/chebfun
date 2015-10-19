% Test file for @chebfun/or.m.

function pass = test_or(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Check the empty case.
f = chebfun(@sin);
pass(1) = isempty(f | chebfun()) && isempty(chebfun() | f);

% Check a few simple examples.
f = chebfun(@(x) sin(x), [-1 -0.5 0.5 1], pref);

g = chebfun(@(x) 0*x);
h = f | g;
ind = find(h.pointValues == 0);
pass(2) = all(feval(h, x) == 1) && (numel(ind) == 1) && ...
    (abs(h.domain(ind)) < 10*vscale(h)*eps);

g = chebfun(@(x) exp(x));
h = f | g;
pass(3) = isequal(h.pointValues, [1 1].') && all(feval(h, x) == 1);

f = chebfun(@(x) [sin(x) cos(x)], pref);
g = chebfun(@(x) [0*x exp(x)], pref);

h = f | g;
h_exact = @(x) [(sin(x) | 0*x) (0*x + 1)];
err = feval(h, x) - h_exact(x);
pass(4) = (norm(err(:), inf) == 0) && isequal(h.pointValues, [1 1 ; 0 1 ; 1 1]);

h = f.' | g.';
h_exact = @(x) [(sin(x) | 0*x) (0*x + 1)].';
err = feval(h, x) - h_exact(x);
pass(5) = (norm(err(:), inf) == 0) && isequal(h.pointValues, [1 1 ; 0 1 ; 1 1]);

% Check error conditions.
try
    f = chebfun(@(x) sin(x));
    g = chebfun(@(x) [sin(x) exp(x)]);
    h = f | g;
    pass(6) = false;
catch ME
    pass(6) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:or:dims');
end

try
    f = chebfun(@(x) sin(x));
    g = chebfun(@(x) cos(x), [-2 7]);
    h = f | g;
    pass(7) = false;
catch ME
    pass(7) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:or:doms');
end

%% Test for singular function:

% define the domain:
dom = [-2 7];

op = @(x) sin(30*x)./((x-dom(1)).*(x-dom(2)));
f = chebfun(op, dom, 'exps', [-1 -1], 'splitting', 'on');
g = chebfun(@(x) 0*x, dom, 'splitting', 'on');
h1 = f | f;
h2 = f | g;

% check values:

% Generate a few random points to use as test values:
x = diff(dom) * rand(100, 1) + dom(1);

fval = feval(h1, x);
err = fval - 1;
pass(8) = ~any( err );

r = roots(f);
pass(9) = ~any( h1(r) );

fval = feval(h2, x);
err = fval - 1;
pass(10) = ~any( err );

r = roots(f);
pass(11) = ~any( h2(r) );

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
g = chebfun(@(x) 0*x, dom([1 3]));
h = logical(f);
hVals = feval(h, x);
hExact = opExact(x);
err = hVals - hExact;
pass(12) = ( ~any(err) ) && all(h.pointValues == [1 1 0].');

end

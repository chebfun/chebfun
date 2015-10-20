% Test file for @chebfun/and.m.

function pass = test_and(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Check the empty case.
f = chebfun(@sin);
pass(1) = isempty(f & chebfun()) && isempty(chebfun() & f);

% Check a few simple examples.
f = chebfun(@(x) sin(x), [-1 -0.5 0.5 1], pref);

g = chebfun(@(x) 0*x);
h = f & g;
pass(2) = all(h.pointValues == 0) && all(feval(h, x) == 0);

g = chebfun(@(x) exp(x));
h = f & g;
ind = find(h.pointValues == 0);
pass(3) = all(feval(h, x) == 1) && (numel(ind) == 1) && ...
    (abs(h.domain(ind)) < 10*vscale(h)*eps);

f = chebfun(@(x) [sin(x) cos(x)], pref);
g = chebfun(@(x) [exp(x) 0*x], pref);

h = f & g;
h_exact = @(x) [(sin(x) & exp(x)) 0*x];
err = feval(h, x) - h_exact(x);
pass(4) = (norm(err(:), inf) == 0) && isequal(h.pointValues, [1 0 ; 0 0 ; 1 0]);

h = f.' & g.';
h_exact = @(x) [(sin(x) & exp(x)) 0*x].';
err = feval(h, x) - h_exact(x);
pass(5) = (norm(err(:), inf) == 0) && isequal(h.pointValues, [1 0 ; 0 0 ; 1 0]);

% Check error conditions.
try
    f = chebfun(@(x) sin(x));
    g = chebfun(@(x) [sin(x) exp(x)]);
    h = f & g;
    pass(6) = false;
catch ME
    pass(6) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:and:dims');
end

try
    f = chebfun(@(x) sin(x));
    g = chebfun(@(x) cos(x), [-2 7]);
    h = f & g;
    pass(7) = false;
catch ME
    pass(7) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:and:doms');
end

%% Test for singular functions:

% define the domain:
dom = [-2 7];

op = @(x) sin(30*x)./((x-dom(1)).*(x-dom(2)));
f = chebfun(op, dom, 'exps', [-1 -1], 'splitting', 'on');
g = chebfun(@(x) x.^2+1, dom, 'splitting', 'on');
h1 = f & f;
h2 = f & g;

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
pass(11) = ~any( h2(r) );

%% Test for function defined on unbounded domain:

% Set the domain:
dom = [1 Inf];
domCheck = [1 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) 1./x;
f = chebfun(op, dom);
g = chebfun(@(x) 0*x, dom);
h = f & g;
hVals = feval(h, x);
pass(12) = ~any( hVals );

end

% Test file for @chebfun/any.m.

function pass = test_any(pref)

if ( nargin <  1 )
    pref = chebfunpref();
end

% Enable breakpoint detection.
pref.splitting = 1;

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;
hvsde = @(x) .5*(sign(x) + 1);

% Check empty case.
pass(1) = ~any(chebfun());
pass(2) = isempty(any(chebfun(), 2));

% Check behavior for any() along continuous dimension (vector output).
f = chebfun(@(x) [sin(x) 0*x exp(x)], [-1 -0.5 0 0.5 1], pref);
pass(2) = isequal(any(f), [1 0 1]);
pass(3) = isequal(any(f.', 2), [1 0 1].');

f.pointValues(3,2) = NaN;
pass(4) = isequal(any(f, 1), [1 0 1]);

f = chebfun(@(x) [0*x hvsde(x) exp(2*pi*1i*x)], [-1 0 1], pref);
pass(5) = isequal(any(f), [0 1 1]);
pass(6) = isequal(any(f.', 2), [0 1 1].');

% Check behavior for any() along discrete dimension (chebfun output).
f = chebfun(@(x) [sin(x) 0*x exp(x)], pref);
g = any(f, 2);
pass(7) = ~g.isTransposed && (numel(g.funs) == 1) && all(feval(g, x) == 1);
g = any(f.', 1);
pass(8) = g.isTransposed && (numel(g.funs) == 1) && all(feval(g, x) == 1);

f = chebfun(@(x) [sin(x) 0*x], pref);
g = any(f, 2);
ind = find(g.pointValues == 0);
pass(9) = ~g.isTransposed && (abs(g.domain(ind)) < 10*vscale(g)*epslevel(g)) ...
    && isequal(g.pointValues, [1 0 1].') && all(feval(g, x) == 1);
g = any(f.', 1);
ind = find(g.pointValues == 0);
pass(10) = g.isTransposed && (abs(g.domain(ind)) < 10*vscale(g)*epslevel(g)) ...
    && isequal(g.pointValues, [1 0 1].') && all(feval(g, x) == 1);

f = chebfun(@(x) [hvsde(x) sin(x).*hvsde(x)], [-1 0 1], pref);
g = any(f, 2);
g_exact = @(x) any([hvsde(x) sin(x).*hvsde(x)], 2);
pass(11) = ~g.isTransposed && isequal(g.pointValues, [0 1 1].') && ...
    all(feval(g, x) == g_exact(x));

% Check error conditions.
try
    any(f, 3)
    pass(12) = false;
catch ME
    pass(12) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:any:dim');
end

%% Test for singular function:

% define the domain:
dom = [-2 7];

op = @(x) sin(30*x)./((x-dom(1)).*(x-dom(2)));
f = chebfun(op, dom, 'exps', [-1 -1], 'splitting', 'on');
h1 = any(f, 1);
h2 = any(f, 2);

% check values:

% Generate a few random points to use as test values:
x = diff(dom) * rand(100, 1) + dom(1);

pass(13) = ~( h1 - 1 );

fval = feval(h2, x);
err = fval - 1;
pass(14) = ~any( err );

%% Test for function defined on unbounded domain:

% Functions on [a inf]:

% Set the domain:
dom = [1 Inf];

op = @(x) 1./x.^2;
f = chebfun(op, dom, 'exps', [0 -2]);
pass(15) = any(f);

g = chebfun(@(x) 0*x, dom);
pass(16) = ~any(g);

end

% Test file for @chebfun/le.m.

function pass = test_le(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;
hvsde = @(x) .5*(sign(x) + 1);

% Check the empty case.
f = chebfun(@(x) sin(x), [-1 -0.5 0 0.5 1], pref);
g = chebfun();
pass(1) = isempty(f <= g)  && isempty(g <= f);

% Check a few simple examples.
g = chebfun(@(x) 0*x + sqrt(2)/2, pref);
h = f <= g;
ind = find(abs(h.domain - pi/4) < 10*vscale(h)*epslevel(h));
pass(2) = ~isempty(ind) && (h.pointValues(ind) == 1) && ...
    all(feval(h, x) - hvsde(sqrt(2)/2 - sin(x)) == 0);

f = chebfun(@(x) exp(x), pref);
g_op = @(x) (exp(0.5) - exp(-0.5))*(x + 0.5) + exp(-0.5);
g = chebfun(g_op, pref);
h = f <= g;
ind1 = find(abs(h.domain - 0.5) < 10*vscale(h)*epslevel(h));
ind2 = find(abs(h.domain + 0.5) < 10*vscale(h)*epslevel(h));

pass(3) = ~isempty(ind1) && (h.pointValues(ind1) == 1) && ...
    ~isempty(ind2) && (h.pointValues(ind2) == 1) && ...
    all(feval(h, x) - hvsde(g_op(x) - exp(x)) == 0);

h = f <= f;
pass(4) = (numel(h.funs) == 1) && all(feval(h, x) == 1);

h = f <= -f;
pass(5) = (numel(h.funs) == 1) && all(feval(h, x) == 0);

%% Check an example where le() does not need to introduce any new breakpoints.
f = chebfun(@(x) exp(x), [-1 -0.5 0 0.5 1], pref);
h = f <= g;
ind1 = find(abs(h.domain - 0.5) < 10*vscale(h)*epslevel(h));
ind2 = find(abs(h.domain + 0.5) < 10*vscale(h)*epslevel(h));

pass(6) = ~isempty(ind1) && (h.pointValues(ind1) == 1) && ...
    ~isempty(ind2) && (h.pointValues(ind2) == 1) && ...
    all(feval(h, x) - hvsde(g_op(x) - exp(x)) == 0);

% Check error conditions.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) exp(2*pi*1i*x), [-1 1], pref);

try
    h = f <= g
    pass(7) = false;
catch ME
    pass(7) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:le:array');
end

try
    h = g <= f
    pass(8) = false;
catch ME
    pass(8) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:le:array');
end

%% Test for singular function:

f = chebfun(@(x) -sin(x)./(x+1), 'exps', [-1 0]);
g = chebfun(@(x) x*0);
h = ( f <= g );
h_vals = feval(h, x);
h_exact = hvsde(x);
pass(9) = ( h.pointValues(1) == 0 ) && all( abs( h_vals - h_exact ) == 0 );

%% Test for function defined on unbounded domain:

% Functions on [a inf]:

% Set the domain:
dom = [-1 Inf];
domCheck = [-1 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

opf = @(x) exp(-x)-1;
opg = @(x) x.*exp(-x);
oph = @(x) opf(x) <= opg(x);

f = chebfun(opf, dom);
g = chebfun(opg, dom);
h = ( f <= g );

hVals = feval(h, x);
hExact = oph(x);
err = hVals - hExact;
pass(10) = ( ~any(err) ) && all(h.pointValues == [0 1 1].');

end

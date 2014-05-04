% Test file for sing/compose.m

function pass = test_compose(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
d = 2;
x = 2*(1-10^(-d)) * rand(100, 1) - (1-10^(-d));

% Composition of two SINGFUNs - OP(F, G):
f = singfun(@(x) sin(x)./(x+1), [-1 0], [], [], [], []);
g = singfun(@(x) cos(x)./(x+1), [-1 0], [], [], [], []);
pref = chebfunpref();
pref.singPrefs.exponents = [-1 0];
h = compose(f, @plus, g, pref);
op = @(x) (sin(x)+cos(x))./(x+1);
hVals = feval(h, x);
hExact = op(x);
err = hVals - hExact;
pass(1) = norm(err, inf) < 1e3*get(h, 'epslevel')*get(h, 'vscale');

% Composition of an operator and a SINGFUN - OP(F)
f = singfun(@(x) sqrt(x+1), [0.5 0], [], [], [], []);
h = compose(f, @sin);
op = @(x) sin(sqrt(x+1));
hVals = feval(h, x);
hExact = op(x);
err = hVals - hExact;
pass(2) = norm(err, inf) < get(h, 'epslevel')*get(h, 'vscale');

% Composition of a SMOOTHFUN and a SINGFUN - G(F)
f = singfun(@(x) sqrt(x+1), [0.5 0], [], [], [], []);
g = chebtech2(@(x) cos(x));
h = compose(f, g);
op = @(x) cos(sqrt(x+1));
hVals = feval(h, x);
hExact = op(x);
err = hVals - hExact;
pass(3) = norm(err, inf) < 1e3*get(h, 'epslevel')*get(h, 'vscale');

end

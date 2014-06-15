% Test file for fractional calculas.

function pass = test_fracCalc(pref)

% Obtain preferences.
if ( nargin == 0 )
    pref = chebfunpref();
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

%% Differetiation:
f = chebfun(@(x)sin(x), [-1 1]);
g = diff(f, 0.3);
h = diff(g, 0.7);
v = diff(f);
pass(1) = norm(h-v, inf) < 1e3*epslevel(f)*hscale(f);

%% Differetiation of a piecewise CHEBFUN:
f = chebfun({@(x)sin(x) @(x)sin(x)}, [-1 0 1]);
g = diff(f, 0.3);
h = diff(g, 0.7);
v = diff(f);
pass(2) = norm(h-v, inf) < 1e4*epslevel(f)*hscale(f);

%% Integration:
f = chebfun(@(x)sin(x), [-1 1]);
g = cumsum(f, 0.3);
h = cumsum(g, 0.7);
v = cumsum(f);
pass(3) = norm(h-v, inf) < 1e1*epslevel(f)*hscale(f);

%% Integration of a piecewise smooth CHEBFUN:
f = chebfun({@(x)sin(x) @(x)sin(x)}, [-1 0 1]);
g = cumsum(f, 0.3);
h = cumsum(g, 0.7);
v = cumsum(f);
pass(4) = norm(h-v, inf) < 1e3*epslevel(f)*hscale(f);

%% Differentiation of a singular function:
f = chebfun(@(x)sin(x).*((1+x).^0.3), [-1 1], 'exps', [0.3 0]);
g = diff(f, 0.3);
h = diff(g, 0.7);
v = diff(f);
err = feval(h-v, xr);
pass(5) = norm(err, inf) < 1e5*epslevel(f)*hscale(f);

%% Differentiation of a singular function represented by a piecewise smooth CHEBFUN:
f = chebfun({@(x)sin(x).*((1+x).^0.3) @(x)sin(x).*((1+x).^0.3)}, [-1 0 1], 'exps', [0.3 0 0 0]);
g = diff(f, 0.3);
h = diff(g, 0.7);
v = diff(f);
err = feval(h-v, xr);
pass(6) = norm(err, inf) < 1e4*epslevel(f)*hscale(f);

%% Integration of of a singular function:
f = chebfun(@(x)sin(x).*((1+x).^0.3), [-1 1], 'exps', [0.3 0]);
g = cumsum(f, 0.3);
h = cumsum(g, 0.7);
v = cumsum(f);
err = feval(h-v, xr);
pass(7) = norm(err, inf) < 1e1*epslevel(f)*hscale(f);

%% Integration of a singular function represented by a piecewise smooth CHEBFUN:
f = chebfun({@(x)sin(x).*((1+x).^0.3) @(x)sin(x).*((1+x).^0.3)}, [-1 0 1], 'exps', [0.3 0 0 0]);
g = cumsum(f, 0.3);
h = cumsum(g, 0.7);
v = cumsum(f);
err = feval(h-v, xr);
pass(8) = norm(err, inf) < 1e1*epslevel(f)*hscale(f);

%% Caputo definition:

f = chebfun(@(x)exp(x), [-2 7]);
g = diff(f, 0.3, 'caputo');
h = diff(g, 0.7, 'caputo');
v = diff(f);
pass(9) = norm(h-v, inf) < 1e3*epslevel(f)*hscale(f);



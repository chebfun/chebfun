% Test file for unbndfun/cumsum.m.

function pass = test_cumsum(pref)

if ( nargin == 1 )
    pref = chebfunpref();
end

% Seed for random number:
seedRNG(6178);

%% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];
domCheck = [-1e2 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) exp(-x.^2);
f = unbndfun(op, dom);
g = cumsum(f);

gVals = feval(g, x);
opg = @(x) sqrt(pi)*erf(x)/2 + sqrt(pi)/2;
gExact = opg(x);
errg = norm(gVals - gExact, inf);
tol = 4e3*get(g,'epslevel').*get(g,'vscale');
pass(1) = errg < tol;

% [TODO]: Revive when log is ready.
% Blow-up function:
% op = @(x) x.^2.*(1-exp(-x.^2));
% pref.singPrefs.exponents = [2 2];
% f = unbndfun(op, dom, [], [], pref);
% g = cumsum(f);
% 
% opg = @(x) x.*exp(-x.^2)/2 + x.^3/3 - sqrt(pi)*erf(x)/4;
% gVals = feval(g, x);
% gExact = opg(x);
% err = gVals - gExact;
% pass(3) = norm(err, inf) < 1e4*get(g,'epslevel').*get(g,'vscale');
pass(2) = 1;

%% Functions on [a inf]:

% Set the domain:
dom = [1 Inf];
domCheck = [1 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) x.*exp(-x);
f = unbndfun(op, dom);
g = cumsum(f);
gVals = feval(g, x);

opg = @(x) -exp(-x).*(x + 1) + 2*exp(-1);
gExact = opg(x);
err = gVals - gExact;
pass(3) = norm(err, inf) < 1e5*get(g,'epslevel').*get(g,'vscale');

% Blow-up function:
op = @(x) 5*x;
pref.singPrefs.exponents = [0 1];
f = unbndfun(op, dom, [], [], pref);
g = cumsum(f);
gVals = feval(g, x);

opg = @(x) 5*x.^2/2 - 5/2 + get(g, 'lval');
gExact = opg(x);
err = gVals - gExact;
pass(4) = norm(err, inf) < 1e2*get(g,'epslevel').*get(g,'vscale');

%% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];
domCheck = [-1e6 -3*pi];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) exp(x);
f = unbndfun(op, dom);
g = cumsum(f);
gVals = feval(g, x);

opg = @(x) exp(x);
gExact = opg(x);
err = gVals - gExact;
pass(5) = norm(err, inf) < 4e3*get(g,'epslevel').*get(g,'vscale');

%% Array-valued function:
op = @(x) [exp(x) x.*exp(x)];
f = unbndfun(op, dom);
g = cumsum(f);
gVals = feval(g, x);

opg = @(x) [exp(x) exp(x).*(x - 1)];
gExact = opg(x);
err = gVals - gExact;
pass(6) = norm(err, inf) < 2e4*max(get(g,'epslevel').*get(g,'vscale'));

%% Test on cumulative sum over the columns
h = cumsum(f, 2);
hVals = feval(h, x);

oph = @(x) [exp(x) exp(x).*(x+1)];
hExact = oph(x);
err = hVals - hExact;
pass(7) = norm(err, inf) < max(get(h,'epslevel').*get(h,'vscale'));

end
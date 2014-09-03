% Test file for unbndfun/diff.

function pass = test_diff(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

singPref = pref;
singPref.blowup = true;

% Seed for random number:
seedRNG(6178);

%% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];
domCheck = [-1e2 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) exp(-x.^2);
f = unbndfun(op, struct('domain', dom));
g = diff(f, 2);
op_g = @(x) 4*x.^2.*exp(-x.^2) - 2*exp(-x.^2);
gVals = feval(g, x);
gExact = op_g(x);
err = gVals - gExact;
pass(1) = norm(err, inf) < 5e1*get(g,'epslevel')*get(g,'vscale');

op = @(x) x.^2.*exp(-x.^2);
f = unbndfun(op, struct('domain', dom));
g = diff(f);
op_g = @(x) 2*x.*exp(-x.^2) - 2*x.^3.*exp(-x.^2);
gVals = feval(g, x);
gExact = op_g(x);
err = gVals - gExact;
pass(2) = norm(err, inf) < 2e1*get(g,'epslevel')*get(g,'vscale');

op = @(x) (1-exp(-x.^2))./x;
f = unbndfun(op, struct('domain', dom));
g = diff(f);
op_g = @(x) 2*exp(-x.^2) + (exp(-x.^2) - 1)./x.^2;
gVals = feval(g, x);
gExact = op_g(x);
err = norm(gVals - gExact,inf);
pass(3) = err < 1e1*get(g,'epslevel')*get(g,'vscale');

op = @(x) x.^2.*(1-exp(-x.^2));
f = unbndfun(op, struct('domain', dom, 'exponents', [2 2]), singPref);
g = diff(f);
op_g = @(x) 2*x.^3.*exp(-x.^2) - 2*x.*(exp(-x.^2) - 1);
gVals = feval(g, x);
gExact = op_g(x);
err = norm(gVals - gExact, inf);
tol = 5e5*get(f,'epslevel')*get(f,'vscale');
pass(4) = err < tol;

%% Functions on [a inf]:

% Set the domain:
dom = [1 Inf];
domCheck = [1 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) exp(-x);
f = unbndfun(op, struct('domain', dom));
gVals = feval(f, x);
gExact = op(x);
err = gVals - gExact;
pass(5) = norm(err, inf) < get(f,'epslevel')*get(f,'vscale');

op = @(x) x.*exp(-x);
f = unbndfun(op, struct('domain', dom));
gVals = feval(f, x);
gExact = op(x);
err = gVals - gExact;
pass(6) = norm(err, inf) < get(f,'epslevel')*get(f,'vscale');

op = @(x) (1-exp(-x))./x;
f = unbndfun(op, struct('domain', dom));
gVals = feval(f, x);
gExact = op(x);
err = gVals - gExact;
pass(7) = norm(err, inf) < get(f,'epslevel')*get(f,'vscale');

op = @(x) 1./x;
f = unbndfun(op, struct('domain', dom));
gVals = feval(f, x);
gExact = op(x);
err = gVals - gExact;
pass(8) = norm(err, inf) < get(f,'epslevel')*get(f,'vscale');

op = @(x) x.*(5+exp(-x.^3));
f = unbndfun(op, struct('domain', dom, 'exponents', [0 1]), pref);
gVals = feval(f, x);
gExact = op(x);
err = norm(gVals - gExact, inf);
tol = 10*get(f,'epslevel')*get(f,'vscale');
pass(9) = err < tol;

%% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];
domCheck = [-1e6 -3*pi];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) exp(x);
f = unbndfun(op, struct('domain', dom));
gVals = feval(f, x);
gExact = op(x);
err = gVals - gExact;
pass(10) = norm(err, inf) < get(f,'epslevel')*get(f,'vscale');

op = @(x) x.*exp(x);
f = unbndfun(op, struct('domain', dom));
gVals = feval(f, x);
gExact = op(x);
err = gVals - gExact;
pass(11) = norm(err, inf) < get(f,'epslevel')*get(f,'vscale');

op = @(x) (1-exp(x))./x;
f = unbndfun(op, struct('domain', dom));
gVals = feval(f, x);
gExact = op(x);
err = gVals - gExact;
pass(12) = norm(err, inf) < get(f,'epslevel')*get(f,'vscale');

op = @(x) 1./x;
f = unbndfun(op, struct('domain', dom));
gVals = feval(f, x);
gExact = op(x);
err = gVals - gExact;
pass(13) = norm(err, inf) < get(f,'epslevel')*get(f,'vscale');

op = @(x) x.*(5+exp(x.^3))./(dom(2)-x);
f = unbndfun(op, struct('domain', dom, 'exponents', [0 -1]), pref);
gVals = feval(f, x);
gExact = op(x);
err = gVals - gExact;
pass(14) = norm(err, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');
    
end

% Test file for unbndfun/changeMap.m.

function pass = test_changeMap(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebfunpref();
end

singPref = pref;
singPref.blowup = true;

% Seed for random number:
seedRNG(6178);

%% Functions on [a inf]:

% Set the domain:
dom = [1 Inf];
domCheck = [1 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% New domain:
domNew = [7 Inf];
xNew = domNew(1) + x - dom(1);

% Exponentially decaying function:
op = @(x) (1-exp(-x))./x;
f = unbndfun(op, struct('domain', dom));
g = changeMap(f, domNew);
fVals = feval(f, x);
gVals = feval(g, xNew);
err = fVals - gVals;
pass(1) = norm(err, inf) < get(f,'epslevel')*get(f,'vscale');

% Blow-up function:
op = @(x) x.*(5+exp(-x.^3));
f = unbndfun(op, struct('domain', dom, 'exponents', [0 1]), singPref);
g = changeMap(f, domNew);
fVals = feval(f, x);
gVals = feval(g, xNew);
err = fVals - gVals;
pass(2) = norm(err, inf) < get(f,'epslevel')*get(f,'vscale');

%% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];
domCheck = [-1e6 -3*pi];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% New domain:
domNew = [-Inf 200];
xNew = domNew(2) + x - dom(2);

op = @(x) x.*exp(x);
f = unbndfun(op, struct('domain', dom));
fVals = feval(f, x);
g = changeMap(f, domNew);
gVals = feval(g, xNew);
err = fVals - gVals;
pass(3) = norm(err, inf) < get(f,'epslevel')*get(f,'vscale');

% Blow-up function:
op = @(x) x.*(5+exp(x.^3))./(dom(2)-x);
f = unbndfun(op, struct('domain', dom, 'exponents', [0 -1]), singPref);
fVals = feval(f, x);
g = changeMap(f, domNew);
gVals = feval(g, xNew);
err = fVals - gVals;
pass(4) = norm(err, inf) < get(f,'epslevel')*get(f,'vscale');

end

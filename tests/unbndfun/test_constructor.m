% Test file for unbndfun constructor and unbndfun/feval.m.

function pass = test_constructor(pref)

if ( nargin == 1 )
    pref = chebpref();
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
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(1) = norm(err, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');

op = @(x) x.^2.*exp(-x.^2);
f = unbndfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(2) = norm(err, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');

op = @(x) (1-exp(-x.^2))./x;
f = unbndfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(3) = norm(err, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');

% Blow-up function:
op = @(x) x.^2.*(1-exp(-x.^2));
pref.singPrefs.exponents = [2 2];
f = unbndfun(op, dom, [], [], pref); 
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(4) = norm(err, inf) < 1e4*get(f,'epslevel')*get(f,'vscale');

%% Functions on [a inf]:

% Set the domain:
dom = [1 Inf];
domCheck = [1 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) exp(-x);
f = unbndfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(5) = norm(err, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');

op = @(x) x.*exp(-x);
f = unbndfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(6) = norm(err, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');

op = @(x) (1-exp(-x))./x;
f = unbndfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(7) = norm(err, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');

op = @(x) 1./x;
f = unbndfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(8) = norm(err, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');

% Blow-up function:
op = @(x) x.*(5+exp(-x.^3));
pref.singPrefs.exponents = [0 1];
f = unbndfun(op, dom, [], [], pref); 
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(9) = norm(err, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');

%% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];
domCheck = [-1e6 -3*pi];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) exp(x);
f = unbndfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(10) = norm(err, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');

op = @(x) x.*exp(x);
f = unbndfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(11) = norm(err, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');

op = @(x) (1-exp(x))./x;
f = unbndfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(12) = norm(err, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');

op = @(x) 1./x;
f = unbndfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(13) = norm(err, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');

% Blow-up function:
op = @(x) x.*(5+exp(x.^3))./(dom(2)-x);
pref.singPrefs.exponents = [0 -1];
f = unbndfun(op, dom, [], [], pref); 
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(14) = norm(err, inf) < 1e1*get(f,'epslevel')*get(f,'vscale');

% Array-valued function:
op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
f = unbndfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(15) = norm(err, inf) < 1e1*max(get(f,'epslevel').*get(f,'vscale'));

%% MISC:

try
    f = unbndfun(@(x) exp(-x.^2), [0 1]);
    pass(16) = fail;
catch ME
    pass(16) = strcmp(ME.identifier, 'CHEBFUN:UNBNDFUN:BoundedDomain');
end
    
end
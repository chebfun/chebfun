% Test file for unbndfun/restrict.m.

function pass = test_restrict(pref)

if ( nargin == 1 )
    pref = chebpref();
end

% Seed for random number:
seedRNG(6178);

%% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];
dom1 = [-1e2 -2];
dom2 = [-2 7];
dom3 = [7 1e2];

% Generate a few random points to use as test values:
x1 = diff(dom1) * rand(100, 1) + dom1(1);
x2 = diff(dom2) * rand(100, 1) + dom2(1);
x3 = diff(dom3) * rand(100, 1) + dom3(1);

op = @(x) x.^2.*exp(-x.^2);
f = unbndfun(op, dom);
g = restrict(f, [-inf -2 7 inf]);
g1Vals = feval(g{1}, x1);
g2Vals = feval(g{2}, x2);
g3Vals = feval(g{3}, x3);
fExact1 = op(x1);
fExact2 = op(x2);
fExact3 = op(x3);
err1 = g1Vals - fExact1;
err2 = g2Vals - fExact2;
err3 = g3Vals - fExact3;
pass(1) = ( norm(err1, inf) < 1e1*get(g{1},'epslevel')*get(g{1},'vscale') ...
    && norm(err2, inf) < 1e1*get(g{2},'epslevel')*get(g{2},'vscale') ...
    && norm(err3, inf) < 1e1*get(g{3},'epslevel')*get(g{3},'vscale') );



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
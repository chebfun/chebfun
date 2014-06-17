% Test file for unbndfun/restrict.m.

function pass = test_restrict(pref)

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
dom1 = [-1e2 -2];
dom2 = [-2 7];
dom3 = [7 1e2];
domRestrict = [-inf -2 7 inf];

% Generate a few random points to use as test values:
x1 = diff(dom1) * rand(100, 1) + dom1(1);
x2 = diff(dom2) * rand(100, 1) + dom2(1);
x3 = diff(dom3) * rand(100, 1) + dom3(1);

% Exponentially-decaying function:
op = @(x) x.^2.*exp(-x.^2);
f = unbndfun(op, struct('domain', dom));
g = restrict(f, domRestrict);
g1Vals = feval(g{1}, x1);
g2Vals = feval(g{2}, x2);
g3Vals = feval(g{3}, x3);
fExact1 = op(x1);
fExact2 = op(x2);
fExact3 = op(x3);
err1 = norm(g1Vals - fExact1, inf);
err2 = norm(g2Vals - fExact2, inf);
err3 = norm(g3Vals - fExact3, inf);
tol1 = 10*get(g{1},'epslevel')*get(g{1},'vscale');
tol2 = 10*get(g{2},'epslevel')*get(g{2},'vscale');
tol3 = 10*pref.eps;
pass(1) = ( err1 < tol1 && err2 < tol2 && err3 < tol3);

% Blow-up function:
op = @(x) x.^2.*(1-exp(-x.^2));
f = unbndfun(op, struct('domain', dom, 'exponents', [2 2]), singPref);
g = restrict(f, domRestrict);

g1Vals = feval(g{1}, x1);
g2Vals = feval(g{2}, x2);
g3Vals = feval(g{3}, x3);
fExact1 = op(x1);
fExact2 = op(x2);
fExact3 = op(x3);
err1 = norm(g1Vals - fExact1, inf);
err2 = norm(g2Vals - fExact2, inf);
err3 = norm(g3Vals - fExact3, inf);
tol1 = 100*get(g{1},'epslevel')*get(g{1},'vscale');
tol2 = 10*get(g{2},'epslevel')*get(g{2},'vscale');
tol3 = 10*get(g{3},'epslevel')*get(g{3},'vscale');
pass(2) = ( err1 < tol1 && err2 < tol2 && err3 < tol3);

%% Functions on [a inf]:

% Set the domain:
dom = [1 Inf];
dom1 = [1 4];
dom2 = [4 23];
dom3 = [23 1e2];
domRestrict = [1 4 23 inf];

% Generate a few random points to use as test values:
x1 = diff(dom1) * rand(100, 1) + dom1(1);
x2 = diff(dom2) * rand(100, 1) + dom2(1);
x3 = diff(dom3) * rand(100, 1) + dom3(1);

% Exponentially-decaying function:
op = @(x) (1-exp(-x))./x;
f = unbndfun(op, struct('domain', dom));
g = restrict(f, domRestrict);

g1Vals = feval(g{1}, x1);
g2Vals = feval(g{2}, x2);
g3Vals = feval(g{3}, x3);
fExact1 = op(x1);
fExact2 = op(x2);
fExact3 = op(x3);
err1 = g1Vals - fExact1;
err2 = g2Vals - fExact2;
err3 = g3Vals - fExact3;
pass(3) = ( norm(err1, inf) < 2*get(g{1},'epslevel')*get(g{1},'vscale') ...
    && norm(err2, inf) < get(g{2},'epslevel')*get(g{2},'vscale') ...
    && norm(err3, inf) < 2e1*get(g{3},'epslevel')*get(g{3},'vscale') );

% Blow-up function:
op = @(x) x.*(5+exp(-x.^3));
f = unbndfun(op, struct('domain', dom, 'exponents', [0 1]), singPref);
g = restrict(f, domRestrict);

g1Vals = feval(g{1}, x1);
g2Vals = feval(g{2}, x2);
g3Vals = feval(g{3}, x3);
fExact1 = op(x1);
fExact2 = op(x2);
fExact3 = op(x3);
err1 = g1Vals - fExact1;
err2 = g2Vals - fExact2;
err3 = g3Vals - fExact3;
pass(4) = ( norm(err1, inf) < 2*get(g{1},'epslevel')*get(g{1},'vscale') ...
    && norm(err2, inf) < get(g{2},'epslevel')*get(g{2},'vscale') ...
    && norm(err3, inf) < 20*get(g{3},'epslevel')*get(g{3},'vscale') );

%% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];
dom1 = [-3e2 -199];
dom2 = [-199 -66];
dom3 = [-66 -5*pi];
domRestrict = [-Inf -199 -66 -5*pi];

% Generate a few random points to use as test values:
x1 = diff(dom1) * rand(100, 1) + dom1(1);
x2 = diff(dom2) * rand(100, 1) + dom2(1);
x3 = diff(dom3) * rand(100, 1) + dom3(1);

% Exponentially-decaying function:
op = @(x) x.*exp(x);
f = unbndfun(op, struct('domain', dom));
g = restrict(f, domRestrict);

g1Vals = feval(g{1}, x1);
g2Vals = feval(g{2}, x2);
g3Vals = feval(g{3}, x3);
fExact1 = op(x1);
fExact2 = op(x2);
fExact3 = op(x3);
err1 = norm(g1Vals - fExact1, inf);
err2 = norm(g2Vals - fExact2, inf);
err3 = norm(g3Vals - fExact3, inf);
tol1 = max(2*get(g{1},'vscale'), eps);
tol2 = max(2*get(g{2},'vscale'), eps);
tol3 = max(1e1*get(g{3},'epslevel')*get(g{3},'vscale'), eps);
pass(5) = ( err1 < tol1 && err2 < tol2 && err3 < tol3);

% Blow-up function:
op = @(x) x.*(5+exp(x.^3))./(dom(2)-x);
f = unbndfun(op, struct('domain', dom, 'exponents', [0 -1]), singPref);
g = restrict(f, domRestrict);

g1Vals = feval(g{1}, x1);
g2Vals = feval(g{2}, x2);
g3Vals = feval(g{3}, x3);
fExact1 = op(x1);
fExact2 = op(x2);
fExact3 = op(x3);
err1 = g1Vals - fExact1;
err2 = g2Vals - fExact2;
err3 = g3Vals - fExact3;
pass(6) = ( norm(err1, inf) < 1e1*get(g{1},'epslevel')*get(g{1},'vscale') ...
    && norm(err2, inf) < 1e1*get(g{2},'epslevel')*get(g{2},'vscale') ...
    && norm(err3, inf) < 1e2*get(g{3},'epslevel')*get(g{3},'vscale') );

% Array-valued function:
op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
f = unbndfun(op, struct('domain', dom));
g = restrict(f, domRestrict);

g1Vals = feval(g{1}, x1);
g2Vals = feval(g{2}, x2);
g3Vals = feval(g{3}, x3);
fExact1 = op(x1);
fExact2 = op(x2);
fExact3 = op(x3);
err1 = abs(g1Vals - fExact1);
err2 = abs(g2Vals - fExact2);
err3 = abs(g3Vals - fExact3);
bound1 = get(g{1},'epslevel').*get(g{1},'vscale');
bound2 = get(g{2},'epslevel').*get(g{2},'vscale');
bound3 = get(g{3},'epslevel').*get(g{3},'vscale');
pass(7) = ( all( max(err1) < max(bound1, eps) ) ...
    && all( max(err2) < max(bound2, eps) ) ...
    && all( max(err3) < max(bound3, eps) ) );
    
end

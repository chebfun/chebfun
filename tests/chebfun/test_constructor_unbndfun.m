% Test file for chebfun constructor for functions defined in unbounded domain.

function pass = test_constructor_unbndfun(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

%% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];
domCheck = [-1e2 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) exp(-x.^2);
f = chebfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(1) = norm(err, inf) < 1e1*epslevel(f)*vscale(f);

op = @(x) x.^2.*exp(-x.^2);
f = chebfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(2) = norm(err, inf) < 1e1*epslevel(f)*vscale(f);

op = @(x) (1-exp(-x.^2))./x;
f = chebfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(3) = norm(err, inf) < 1e1*epslevel(f)*vscale(f);

% Blow-up function:
op = @(x) x.^2.*(1-exp(-x.^2));
f = chebfun(op, dom, 'exps', [2 2]);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(4) = norm(err, inf) < 1e4*epslevel(f)*vscale(f);

%% Functions on [a inf]:

% Set the domain:
dom = [1 Inf];
domCheck = [1 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) exp(-x);
f = chebfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(5) = norm(err, inf) < 1e1*epslevel(f)*vscale(f);

op = @(x) x.*exp(-x);
f = chebfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(6) = norm(err, inf) < 1e1*epslevel(f)*vscale(f);

op = @(x) (1-exp(-x))./x;
f = chebfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(7) = norm(err, inf) < 1e1*epslevel(f)*vscale(f);

op = @(x) 1./x;
f = chebfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(8) = norm(err, inf) < 1e1*epslevel(f)*vscale(f);

% Blow-up function:
op = @(x) x.*(5+exp(-x.^3));
f = chebfun(op, dom, 'exps', [0 1]); 
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(9) = norm(err, inf) < 1e1*epslevel(f)*vscale(f);

%% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];
domCheck = [-1e6 -3*pi];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) exp(x);
f = chebfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(10) = norm(err, inf) < 1e1*epslevel(f)*vscale(f);

op = @(x) x.*exp(x);
f = chebfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(11) = norm(err, inf) < 1e1*epslevel(f)*vscale(f);

op = @(x) (1-exp(x))./x;
f = chebfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(12) = norm(err, inf) < 1e1*epslevel(f)*vscale(f);

op = @(x) 1./x;
f = chebfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(13) = norm(err, inf) < 1e1*epslevel(f)*vscale(f);

% Blow-up function:
op = @(x) x.*(5+exp(x.^3))./(dom(2)-x);
f = chebfun(op, dom, 'exps', [0 -1]); 
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(14) = norm(err, inf) < 1e1*epslevel(f)*vscale(f);

% Array-valued function:
op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
f = chebfun(op, dom);
fVals = feval(f, x);
fExact = op(x);
err = fVals - fExact;
pass(15) = norm(err, inf) < 1e1*max(epslevel(f).*vscale(f));

end

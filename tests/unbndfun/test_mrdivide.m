% Test file for unbndfun/mrdivide.m

function pass = test_mrdivide(pref)

% Get preferences.
if (nargin < 1)
    pref = chebfunpref();
end

% Set the domain:
dom = [-Inf 3*pi];
domCheck = [-1e6 3*pi];

%% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% Case 1: X*B = A, where B is a scalar and A is an array-valued UNBNDFUN
%         ==> X = A/B, i.e. UNBNDFUN / double

op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
A = unbndfun(op, struct('domain', dom));
B = 3;
X = A/B; 
opExact = @(x) [exp(x)/3 x.*exp(x)/3 (1-exp(x))./(3*x)];
XVals = feval(X, x);
XExact = opExact(x);
err = XVals - XExact;
pass(1) = norm(err, inf) < 2*max(get(X,'epslevel').*get(X,'vscale'));

%% Case 2: X*B = A, where B is a numerical matrix and A is an array-valued 
%         UNBNDFUN ==> X = A/B, i.e. UNBNDFUN / numerical matrix

op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
A = unbndfun(op, struct('domain', dom));
B = rand(3,3);
X = A/B; 
res = X*B - A;
err = feval(res, x);
pass(2) = norm(err(:), inf) < 1e1*max(get(X,'epslevel').*get(X,'vscale'));

end

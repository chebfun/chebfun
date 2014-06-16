% Test file for unbndfun/mldivide.m

function pass = test_mldivide(pref)

% Get preferences.
if (nargin < 1)
    pref = chebfunpref();
end

% Set the domain:
dom = [-Inf 3*pi];
domCheck = [-1e6 3*pi];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% A*X = B, where both A and B are UNBNDFUNs ==> X = A\B

opA = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
A = unbndfun(opA, struct('domain', dom));
opB = @(x) [(2*x+1).*exp(x) exp(x) 2*(1-exp(x))./x];
B = unbndfun(opB, struct('domain', dom));
X = A\B;
res = A*X - B;
err = feval(res, x);
pass(1) = norm(err(:), inf) < max([get(A,'epslevel').*get(A,'vscale') ...
    get(B,'epslevel').*get(B,'vscale')]);

end

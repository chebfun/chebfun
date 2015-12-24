% Test file for @chebfun/fix.m

function pass = test_fix(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Check empty case.
pass(1) = isempty(fix(chebfun()));

% Check "ordinary" cases.
f = chebfun(@(x) 2*sin(x), [-1 -0.5 0.5 1], pref);
g = fix(f);
g_exact = @(x) fix(2*sin(x));
pass(2) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*eps);

f = chebfun(@(x) exp(x), [-1 -0.5 0.5 1], pref);
g = fix(f);
g_exact = @(x) fix(exp(x));
pass(3) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*eps);

% Check complex-valued functions.
f2 = 1i*f;
g = fix(f2);
g_exact = @(x) fix(1i*exp(x));
pass(4) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*eps);

f3 = f + 1i*f;
g = fix(f3);
g_exact = @(x) fix(exp(x) + 1i*exp(x));
pass(5) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*eps);

% Check array-valued function.
f = chebfun(@(x) [2*sin(x) exp(x)], [-1 -0.5 0.5 1], pref);
g = fix(f);
g_exact = @(x) fix([2*sin(x) exp(x)]);
err = feval(g, x) - g_exact(x);
pass(6) = (norm(err(:), inf) <= 10*vscale(g)*eps);

% Check error conditions.
% [TODO]: what to do here
% try
%     f = chebfun(@(x) sin(x), [-1 -0.5 0 0.5 1], pref);
%     f.impulses(3,1,2) = 1;
%     g = fix(f);
%     pass(7) = false;
% catch ME
%     pass(7) = strcmp(ME.identifier, 'CHEBFUN:fix:inf');
% end

% Test for function defined on unbounded domain:

% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];
domCheck = [-1e6 -3*pi];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) x.*exp(x);
f = chebfun(op, dom);
fVals = feval(f, x);
opExact = @(x) fix(op(x));
fExact = opExact(x);
err = fVals - fExact;
pass(7) = norm(err, inf) < eps*vscale(f);

end

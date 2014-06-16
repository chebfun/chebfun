% Test file for @chebfun/floor.m

function pass = test_floor(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Check empty case.
pass(1) = isempty(floor(chebfun()));

% Check "ordinary" cases.
f = chebfun(@(x) sin(x), [-1 -0.5 0.5 1], pref);
g = floor(f);
g_exact = @(x) sign(x)/2 - 1/2;
pass(2) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*epslevel(g));

f = chebfun(@(x) exp(x), [-1 -0.5 0.5 1], pref);
g = floor(f);
g_exact = @(x) floor(exp(x));
pass(3) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*epslevel(g));

% Check complex-valued functions.
f2 = 1i*f;
g = floor(f2);
g_exact = @(x) floor(1i*exp(x));
pass(4) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*epslevel(g));

f3 = f + 1i*f;
g = floor(f3);
g_exact = @(x) floor(exp(x) + 1i*exp(x));
pass(5) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*epslevel(g));

% Check array-valued function.
f = chebfun(@(x) [sin(x) exp(x)], [-1 -0.5 0.5 1], pref);
g = floor(f);
g_exact = @(x) floor([sin(x) exp(x)]);
err = feval(g, x) - g_exact(x);
pass(6) = (norm(err(:), inf) <= 10*vscale(g)*epslevel(g));

% Check error conditions.
% [TODO]: what to do?
% try
%     f = chebfun(@(x) sin(x), [-1 -0.5 0 0.5 1], pref);
%     f.impulses(3,1,2) = 1;
%     g = floor(f);
%     pass(7) = false;
% catch ME
%     pass(7) = strcmp(ME.identifier, 'CHEBFUN:floor:inf');
% end

% Test for function defined on unbounded domain:

% Functions on [a inf]:

% Set the domain:
dom = [1 Inf];
domCheck = [1 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) 1./x.^3;
f = chebfun(op, dom);
g = floor(f);
opExact = @(x) floor(op(x));
fVals = feval(g, x);
fExact = opExact(x);
err = norm(fVals - fExact, inf);
tol = 1e2*epslevel(f)*vscale(f);
pass(7) = err < tol;

end

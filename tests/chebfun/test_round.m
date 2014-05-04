% Test file for @chebfun/round.m

function pass = test_round(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Check empty case.
pass(1) = isempty(round(chebfun()));

% Check "ordinary" cases.
f = chebfun(@(x) sin(x), [-1 -0.5 0.5 1], pref);
g = round(f);
g_exact = @(x) round(sin(x));
pass(2) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*epslevel(g));

f = chebfun(@(x) exp(x), [-1 -0.5 0.5 1], pref);
g = round(f);
g_exact = @(x) round(exp(x));
pass(3) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*epslevel(g));

% Check complex-valued functions.
f2 = 1i*f;
g = round(f2);
g_exact = @(x) round(1i*exp(x));
pass(4) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*epslevel(g));

f3 = f + 1i*f;
g = round(f3);
g_exact = @(x) round(exp(x) + 1i*exp(x));
pass(5) = (norm(feval(g, x) - g_exact(x), inf) <= 10*vscale(g)*epslevel(g));

% Check array-valued function.
f = chebfun(@(x) [sin(x) exp(x)], [-1 -0.5 0.5 1], pref);
g = round(f);
g_exact = @(x) round([sin(x) exp(x)]);
err = feval(g, x) - g_exact(x);
pass(6) = (norm(err(:), inf) <= 10*vscale(g)*epslevel(g));

% Check error conditions.
% [TODO]: What to do here
% try
%     f = chebfun(@(x) sin(x), [-1 -0.5 0 0.5 1], pref);
%     f.impulses(3,1,2) = 1;
%     g = round(f);
%     pass(7) = false;
% catch ME
%     pass(7) = strcmp(ME.identifier, 'CHEBFUN:round:inf');
% end

end

% Test file for singfun/plus.m
% This has been adapted from bndfun test for plus.
function pass = test_plus(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(666);
x = 2 * rand(100, 1) - 1;

% A random number to use as an arbitrary additive constant.
alpha = -0.194758928283640 + 0.075474485412665i;

%%
% Check operation in the case of empty arguments.
f = singfun();
data.exponents = [-1, 0];
g = singfun(@(x) 1./(1+x), data, pref);
pass(1) = (isempty(f + f) && isempty(f + g) && isempty(g + f));

% Addition of smooth SINGFUNs should not return a SINGFUN
f = singfun(@(x) sin(x));
g = singfun(@(x) cos(x));
pass(2) = ~isa(f+g, 'singfun');


% SMOOTHFUN + SINGFUN
f = smoothfun.constructor(@(x) sin(x));
g = singfun(@(x) cos(x));
pass(3) = isa(f + g, 'smoothfun') && isa(g + f, 'smoothfun');
%%
% Check addition with scalars.

fh = @(x) 1./((1+x).*(1-x));
data.exponents = [-1, -1];
f = singfun(fh, data, pref);
pass(4:5) = test_add_function_to_scalar(f, fh, alpha, x);

%%
% Check addition of two singfun objects.

fh = @(x) zeros(size(x));
data.exponents = [];
data.singType = {'none', 'none'};
f = singfun(fh, data, pref);
pass(6:7) = test_add_function_to_function(f, fh, f, fh, x);

fh = @(x) sin(pi*x)./(1-x);
f = singfun(fh, [], pref);

gh = @(x) cos(pi*x)./(1-x);
g = singfun(gh, [], pref);
pass(8:9) = test_add_function_to_function(f, fh, g, gh, x);

gh = @(x) cos(1e2*x);
g = singfun(gh, [], pref);
pass(10:11) = test_add_function_to_function(f, fh, g, gh, x);

gh = @(t) sinh(t*exp(2*pi*1i/6));
g = singfun(gh, [], pref);
pass(12:13) = test_add_function_to_function(f, fh, g, gh, x);

%%
% Check that direct construction and PLUS give comparable results.

f = singfun(@(x) x, [], pref);
g = singfun(@(x) cos(x) - 1, [], pref);
h1 = f + g;
h2 = singfun(@(x) x + cos(x) - 1, [], pref);
tol = 10*max(get(h1, 'vscale')*get(h1, 'epslevel'), ...
    get(h2, 'vscale')*get(h2, 'epslevel'));
pass(14) = norm(feval(h1, x) - feval(h2, x), inf) < tol;

end

% Test the addition of a SINGFUN F, specified by Fh, to a scalar C using
% a grid of points X in [a  b] for testing samples.
function result = test_add_function_to_scalar(f, fh, c, x)
    g1 = f + c;
    g2 = c + f;
    result(1) = isequal(g1, g2);
    g_exact = @(x) fh(x) + c;
    result(2) = norm(feval(g1, x) - g_exact(x), inf) < 2e2*get(f, 'epslevel');
end

% Test the addition of two SINGFUN objects F and G, specified by FH and
% GH, using a grid of points X in [-1  1] for testing samples.
function result = test_add_function_to_function(f, fh, g, gh, x)
    h1 = f + g;
    h2 = g + f;
    result(1) = isequal(h1, h2);
    h_exact = @(x) fh(x) + gh(x);
    result(2) = norm(feval(h1, x) - h_exact(x), inf) <= 1e3*get(f, 'epslevel');
end

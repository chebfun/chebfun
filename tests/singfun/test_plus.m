% Test file for singfun/plus.m
% This has been adapted from bndfun test for plus.
function pass = test_plus(pref)

% Get preferences.
if ( nargin < 1 )
    pref = singfun.pref.singfun;
end

% Generate a few random points to use as test values.
seedRNG(666);
x = 2 * rand(100, 1) + -1;

% A random number to use as an arbitrary additive constant.
alpha = -0.194758928283640 + 0.075474485412665i;

pass = zeros(1, 12); % Pre-allocate pass vector

%%
% Check operation in the case of empty arguments.
f = singfun();
g = singfun(@(x) 1./(1+x), [-1, 0]);
pass(1) = (isempty(f + f) && isempty(f + g) && isempty(g + f));

%%
% Check addition with scalars.

fh = @(x) 1./((1+x).*(1-x));
f = singfun(fh, [-1, -1]);
tol = pref.eps;
pass(2:3) = test_add_function_to_scalar(f, fh, alpha, x, 1e3*tol);

%%
% Check addition of two singfun objects.

fh = @(x) zeros(size(x));
f = singfun(fh, [], {'none', 'none'}, pref);
pass(4:5) = test_add_function_to_function(f, fh, f, fh, x, 1e3*tol);

fh = @(x) sin(pi*x)./(1-x);
f = singfun(fh, [], [], pref);

gh = @(x) cos(pi*x)./(1-x);
g = singfun(gh, [], [], pref);
pass(6:7) = test_add_function_to_function(f, fh, g, gh, x, 1e3*tol);

gh = @(x) cos(1e4*x);
g = singfun(gh, [], [], pref);
pass(8:9) = test_add_function_to_function(f, fh, g, gh, x, 1e3*tol);

gh = @(t) sinh(t*exp(2*pi*1i/6));
g = singfun(gh, [], [], pref);
pass(10:11) = test_add_function_to_function(f, fh, g, gh, x, 1e3*tol);

%%
% Check that direct construction and PLUS give comparable results.

f = singfun(@(x) x, [], [], pref);
g = singfun(@(x) cos(x) - 1, [], [], pref);
h1 = f + g;
h2 = singfun(@(x) x + cos(x) - 1, [], [], pref);
pass(12) = norm(feval(h1, x) - feval(h2, x), inf) < 1e3*tol;


end

% Test the addition of a SINGFUN F, specified by Fh, to a scalar C using
% a grid of points X in [a  b] for testing samples.
function result = test_add_function_to_scalar(f, fh, c, x, tol)
    g1 = f + c;
    g2 = c + f;
    result(1) = isequal(g1, g2);
    g_exact = @(x) fh(x) + c;
    result(2) = norm(feval(g1, x) - g_exact(x), inf) < tol;
end

% Test the addition of two SINGFUN objects F and G, specified by FH and
% GH, using a grid of points X in [-1  1] for testing samples.
function result = test_add_function_to_function(f, fh, g, gh, x, tol)
    h1 = f + g;
    h2 = g + f;
    result(1) = isequal(h1, h2);
    h_exact = @(x) fh(x) + gh(x);
    result(2) = norm(feval(h1, x) - h_exact(x), inf) <= tol;
end

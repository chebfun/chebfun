% Test file for funcheb2/minus.

function pass = test_minus(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb.pref;
end

% Generate a few random points to use as test values.
rngstate = rng();
rng(6178);
x = 2 * rand(100, 1) - 1;

% A random number to use as an arbitrary additive constant.
alpha = randn() + 1i*randn();

%%
% Check operation in the face of empty arguments.

f = funcheb2();
g = funcheb2(@(x) x, pref);
pass(1) = (isempty(f - f) && isempty(f - g) && isempty(g - f));

%%
% Check addition with scalars.

f_op = @(x) sin(x);
f = funcheb2(f_op, pref);
pass(2:3) = test_sub_function_and_scalar(f, f_op, alpha, x);

%%
% Check addition of two funcheb2 objects.

f_op = @(x) zeros(size(x));
f = funcheb2(f_op, pref);
pass(4:5) = test_sub_function_and_function(f, f_op, f, f_op, x);

f_op = @(x) exp(x) - 1;
f = funcheb2(f_op, pref);

g_op = @(x) 1./(1 + x.^2);
g = funcheb2(g_op, pref);
pass(6:7) = test_sub_function_and_function(f, f_op, g, g_op, x);

g_op = @(x) cos(1e4*x);
g = funcheb2(g_op, pref);
pass(8:9) = test_sub_function_and_function(f, f_op, g, g_op, x);

g_op = @(t) sinh(t*exp(2*pi*1i/6));
g = funcheb2(g_op, pref);
pass(10:11) = test_sub_function_and_function(f, f_op, g, g_op, x);

%%
% Check operation for vectorized funcheb2 objects.

f_op = @(x) [zeros(size(x)) zeros(size(x)) zeros(size(x))];
f = funcheb2(f_op, pref);
pass(12:13) = test_sub_function_and_function(f, f_op, f, f_op, x);

f_op = @(x) [sin(x) cos(x) exp(x)];
f = funcheb2(f_op, pref);
pass(14:15) = test_sub_function_and_scalar(f, f_op, alpha, x);

g_op = @(x) [cosh(x) airy(1i*x) sinh(x)];
g = funcheb2(g_op, pref);
pass(16:17) = test_sub_function_and_function(f, f_op, g, g_op, x);

% This should fail with a dimension mismatch error.
g_op = @(x) sin(x);
g = funcheb2(g_op, pref);
try
    h = f - g; %#ok<NASGU>
    pass(18) = false;
catch ME
    pass(18) = strcmp(ME.message, 'Matrix dimensions must agree.');
end

%%
% Check that direct construction and MINUS give comparable results.

tol = 10*eps;
f = funcheb2(@(x) x, pref);
g = funcheb2(@(x) cos(x) - 1, pref);
h1 = f - g;
h2 = funcheb2(@(x) x - (cos(x) - 1), pref);
pass(19) = norm(h1.values - h2.values, 'inf') < tol;

%%
% Restore the RNG state.

rng(rngstate);

end

% Test the subtraction of a FUNCHEB2 F, specified by F_OP, to and from a scalar
% ALPHA using a grid of points X in [-1  1] for testing samples.
function result = test_sub_function_and_scalar(f, f_op, alpha, x)
    g1 = f - alpha;
    g2 = alpha - f;
    result(1) = isequal(g1, -g2);
    g_exact = @(x) f_op(x) - alpha;
    result(2) = norm(feval(g1, x) - g_exact(x), 'inf') < 10*g1.epslevel;
end

% Test the subraction of two FUNCHEB2 objects F and G, specified by F_OP and
% G_OP, using a grid of points X in [-1  1] for testing samples.
function result = test_sub_function_and_function(f, f_op, g, g_op, x)
    h1 = f - g;
    h2 = g - f;
    result(1) = isequal(h1, -h2);
    h_exact = @(x) f_op(x) - g_op(x);
    norm(feval(h1, x) - h_exact(x), 'inf');
    result(2) = norm(feval(h1, x) - h_exact(x), 'inf') < 10*h1.epslevel;
end

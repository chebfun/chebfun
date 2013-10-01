% Test file for @chebfun/mtimes.m.

function pass = test_mtimes(pref)

% Get preferences.
if (nargin < 1)
    pref = chebfun.pref;
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Random numbers to use as arbitrary multiplicative constants.
alpha = -0.194758928283640 + 0.075474485412665i;

% Check operation for empty inputs.
f = chebfun(@sin, [-1 1], pref);
pass(1) = isempty(f*[]);
pass(2) = isempty(f*chebfun());

% Turn on splitting, since we'll need it for the rest of the tests.
pref.chebfun.splitting = 1;

% Test multiplication by scalars.
f1_op = @(x) sin(x).*abs(x - 0.1);
f1 = chebfun(f1_op, pref);
pass(3:4) = test_mult_function_by_scalar(f1, f1_op, alpha, x);

% Test use of mtimes() to compute products.
f = chebfun({@(x) sin(2*pi*x), @(x) cos(2*pi*x)}, [-1 0 1]);
g = chebfun({@(x) exp(2*pi*1i*x), @(x) cos(2*pi*x)}, [-1 0 1]);
tol = 10*max(f.vscale*f.epslevel, g.vscale*g.epslevel);
pass(5) = abs(g.'*f - (0.5 + 0.5i)) < tol;

% Test operation for array-valued chebfuns.
f2_op = @(x) [sin(x).*abs(x - 0.1)  exp(x)];
f2 = chebfun(f2_op, pref);
pass(6:7) = test_mult_function_by_scalar(f2, f2_op, alpha, x);

f = chebfun({@(x) [sin(2*pi*x) sin(2*pi*x)], ...
    @(x) [cos(2*pi*x) cos(4*pi*x)]}, [-1 0 1]);
g = chebfun({@(x) [sin(4*pi*x) sin(4*pi*x)], ...
    @(x) [cos(2*pi*x) cos(4*pi*x)]}, [-1 0 1]);
err = g.'*f - [0.5 0; 0 0.5];
tol = 10*max(f.vscale*f.epslevel, g.vscale*g.epslevel);
pass(8) = norm(err(:), inf) < tol;

A = randn(2, 2);
g = f2*A;
g_exact = @(x) [sin(x).*abs(x - 0.1)  exp(x)]*A;
err = abs(feval(g, x) - g_exact(x));
pass(9) = max(err(:)) < 10*vscale(g)*epslevel(g);

% Test error conditions.
try
    h = f*'X';
    pass(10) = false;
catch ME
    pass(10) = strcmp(ME.identifier, 'CHEBFUN:mtimes:unknown');
end

try
    h = f*f1;
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.identifier, 'CHEBFUN:mtimes:dims');
end

try
    h = f*g;
    pass(12) = false;
catch ME
    pass(12) = strcmp(ME.identifier, 'CHEBFUN:mtimes:dims');
end

try
    h = f*g.';
    pass(13) = false;
catch ME
    pass(13) = strcmp(ME.identifier, 'CHEBFUN:mtimes:colTimesRow');
end

end

% Test the multiplication of a chebfun F, specified by F_OP, by a scalar ALPHA
% using a grid of points X in [-1  1] for testing samples.
function result = test_mult_function_by_scalar(f, f_op, alpha, x)
    g1 = f * alpha;
    g2 = alpha * f;
    result(1) = isequal(g1, g2);
    g_exact = @(x) f_op(x) * alpha;
    err = feval(g1, x) - g_exact(x);
    result(2) = norm(err(:), inf) < 10*max(g1.vscale.*g1.epslevel);
end


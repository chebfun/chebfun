% Test file for @chebfun/minus.m.

function pass = test_minus(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% A random number to use as an arbitrary additive constant.
alpha = -0.194758928283640 + 0.075474485412665i;

% Check behavior for empty arguments.
f = chebfun(@(x) sin(x), pref);
g = chebfun();
pass(1) = isempty(f - []);
pass(2) = isempty(f - g);

% Turn on splitting, since we'll need it for the rest of the tests.
pref.splitting = 1;

% Test addition with scalars.
f1_op = @(x) sin(x).*abs(x - 0.1);
f1 = chebfun(f1_op, pref);
%keyboard;
pass(3:4) = test_sub_function_and_scalar(f1, f1_op, alpha, x);

% Test subition of two chebfun objects.
g1_op = @(x) cos(x).*sign(x + 0.2);
g1 = chebfun(g1_op, pref);
pass(5:6) = test_sub_function_and_function(f1, f1_op, g1, g1_op, x);

% Test operation for array-valued chebfuns.
f2_op = @(x) [sin(x).*abs(x - 0.1)  exp(x)];
f2 = chebfun(f2_op, pref);
pass(7:8) = test_sub_function_and_scalar(f2, f2_op, alpha, x);

g2_op = @(x) [cos(x).*sign(x + 0.2) tan(x)];
g2 = chebfun(g2_op, pref);
pass(9:10) = test_sub_function_and_function(f2, f2_op, g2, g2_op, x);

% Test operation for transposed chebfuns.
pass(11:12) = test_sub_function_and_scalar(f1.', @(x) f1_op(x).', alpha, x);
pass(13:14) = test_sub_function_and_function(f1.', @(x) f1_op(x).', ...
    g1.', @(x) g1_op(x).', x);

% Check error conditions.
try
    h = f1 - uint8(128);
    pass(15) = strcmp(ME.identifier, 'CHEBFUN:plus:unknown')
catch ME
    pass(15) = true;
end

try
    h = f1 - g1.'
    pass(16) = strcmp(ME.identifier, 'CHEBFUN:plus:matdim')
catch ME
    pass(16) = true;
end

%% Tests for function defined on unbounded domain:

% Functions on [-inf b]:

% Set the domain:
dom = [-Inf 3*pi];
domCheck = [-1e6 3*pi];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

opf = @(x) x.*exp(x);
opg = @(x) (1-exp(x))./x;
oph = @(x) x.*exp(x) - (1-exp(x))./x;
f = chebfun(opf, dom);
g = chebfun(opg, dom);
h = f - g;

hVals = feval(h, x);
hExact = oph(x);
err = hVals - hExact;
pass(17) = norm(err, inf) < 1e1*eps*get(h,'vscale');


end

% Test the subtraction of a CHEBFUN F, specified by F_OP, to and from a scalar
% ALPHA using a grid of points X in the domain of F for testing samples.
function result = test_sub_function_and_scalar(f, f_op, alpha, x)
    g1 = f - alpha;
    g2 = alpha - f;
    result(1) = isequal(g1, -g2);
    g_exact = @(x) f_op(x) - alpha;
    result(2) = norm(feval(g1, x) - g_exact(x), inf) < 100*vscale(g1)*eps;
end

% Test the subraction of two CHEBFUN objects F and G, specified by F_OP and
% G_OP, using a grid of points X in the domains of F and G for testing samples.
function result = test_sub_function_and_function(f, f_op, g, g_op, x)
    h1 = f - g;
    h2 = g - f;
    result(1) = isequal(h1, -h2);
    h_exact = @(x) f_op(x) - g_op(x);
    norm(feval(h1, x) - h_exact(x), inf);
    result(2) = norm(feval(h1, x) - h_exact(x), inf) < 100*vscale(h1)*eps;
end

% Test file for @chebfun/plus.m.

function pass = test_plus(pref)

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
pass(1) = isempty(f + []);
pass(2) = isempty(f + g);

% Turn on splitting, since we'll need it for the rest of the tests.
pref.enableBreakpointDetection = 1;

%% Test addition with scalars.
f1_op = @(x) sin(x).*abs(x - 0.1);
f1 = chebfun(f1_op, pref);
pass(3:4) = test_add_function_to_scalar(f1, f1_op, alpha, x);

%% Test addition of two chebfun objects.
g1_op = @(x) cos(x).*sign(x + 0.2);
g1 = chebfun(g1_op, pref);
pass(5:6) = test_add_function_to_function(f1, f1_op, g1, g1_op, x);

% Test operation for array-valued chebfuns.
f2_op = @(x) [sin(x).*abs(x - 0.1)  exp(x)];
f2 = chebfun(f2_op, pref);
pass(7:8) = test_add_function_to_scalar(f2, f2_op, alpha, x);

g2_op = @(x) [cos(x).*sign(x + 0.2) tan(x)];
g2 = chebfun(g2_op, pref);
pass(9:10) = test_add_function_to_function(f2, f2_op, g2, g2_op, x);

% Test operation for transposed chebfuns.
pass(11:12) = test_add_function_to_scalar(f1.', @(x) f1_op(x).', alpha, x);
pass(13:14) = test_add_function_to_function(f1.', @(x) f1_op(x).', ...
    g1.', @(x) g1_op(x).', x);

% Check error conditions.
try
    h = f1 + uint8(128);
    pass(15) = strcmp(ME.identifier, 'CHEBFUN:plus:unknown');
catch ME
    pass(15) = true;
end

try
    h = f1 + g1.';
    pass(16) = strcmp(ME.identifier, 'CHEBFUN:plus:matdim');
catch ME
    pass(16) = true;
end

% Test addition of array-valued scalar to array-valued chebfun.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], pref);
g = f + [1 2 3];
g_exact = @(x) [(1 + sin(x)) (2 + cos(x)) (3 + exp(x))];
err = feval(g, x) - g_exact(x);
pass(17) = norm(err(:), inf) < 10*max(g.vscale*g.epslevel);

% Test scalar expansion in chebfun argument.
f = chebfun(@(x) sin(x), pref);
g = f + [1 2 3];
g_exact = @(x) [(1 + sin(x)) (2 + sin(x)) (3 + sin(x))];
err = feval(g, x) - g_exact(x);
pass(18) = isequal(size(g, 2), 3) && norm(err(:), inf) < ...
    10*max(g.vscale*g.epslevel);

%% QUASIMATRIX

% Test operation for quasimatrices.
f2_op = @(x) [sin(x).*abs(x - 0.1)  exp(x)];
f2 = chebfun(f2_op, pref);
f2q = quasimatrix(f2_op, pref);
pass(19:20) = test_add_function_to_scalar(f2q, f2_op, alpha, x);

g2_op = @(x) [cos(x).*sign(x + 0.2) tan(x)];
g2 = chebfun(g2_op, pref);
g2q = quasimatrix(g2_op, pref);
pass(21:22) = test_add_function_to_function(f2q, f2_op, g2q, g2_op, x);
pass(23:24) = test_add_function_to_function(f2, f2_op, g2q, g2_op, x);
pass(25:26) = test_add_function_to_function(f2q, f2_op, g2, g2_op, x);

% Test addition of array-valued scalar to quasimatrices.
f = quasimatrix(@(x) [sin(x) cos(x) exp(x)], pref);
g = f + [1 2 3];
g_exact = @(x) [(1 + sin(x)) (2 + cos(x)) (3 + exp(x))];
err = feval(g, x) - g_exact(x);
pass(27) = norm(err(:), inf) < 10*max(vscale(g)*epslevel(g));

%% Test on singular function: piecewise smooth chebfun - splitting on.

dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

pow = -1;
op1 = @(x) (x - dom(2)).^pow.*sin(100*x);
op2 = @(x) (x - dom(2)).^pow.*cos(300*x);
pref.singPrefs.exponents = [0 pow];
pref.enableBreakpointDetection = 1;
f = chebfun(op1, dom, pref);
g = chebfun(op2, dom, pref);
h = f + g;
vals_h = feval(h, x);
op = @(x)  (x - dom(2)).^pow.*(sin(100*x)+cos(300*x));
h_exact = op(x);
pass(28) = ( norm(vals_h-h_exact, inf) < max(get(f, 'epslevel'), get(g, 'epslevel'))*...
    norm(h_exact, inf) );

%% Test for function defined on unbounded domain:

% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];
domCheck = [-1e2 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

opf = @(x) exp(-x.^2);
opg = @(x) x.^2.*exp(-x.^2);
oph = @(x) exp(-x.^2) + x.^2.*exp(-x.^2);
f = chebfun(opf, dom);
g = chebfun(opg, dom);
h = f + g;
hVals = feval(h, x);
hExact = oph(x);
err = hVals - hExact;
pass(29) = norm(err, inf) < get(h,'epslevel').*get(h,'vscale');

end

% Test the addition of a chebfun F, specified by F_OP, to a scalar ALPHA using
% a grid of points X in the domain of F for testing samples.
function result = test_add_function_to_scalar(f, f_op, alpha, x)
    g1 = f + alpha;
    g2 = alpha + f;
    result(1) = isequal(g1, g2);
    g_exact = @(x) f_op(x) + alpha;
    result(2) = norm(feval(g1, x) - g_exact(x), inf) < 10*vscale(g1)*epslevel(g1);
end

% Test the addition of two chebfun objects F and G, specified by F_OP and
% G_OP, using a grid of points X in the domain of F and G for testing samples.
function result = test_add_function_to_function(f, f_op, g, g_op, x)
    h1 = f + g;
    h2 = g + f;
    result(1) = isequal(h1, h2);
    h_exact = @(x) f_op(x) + g_op(x);
    norm(feval(h1, x) - h_exact(x), inf);
    result(2) = norm(feval(h1, x) - h_exact(x), inf) < 10*vscale(h1)*epslevel(h1);
end

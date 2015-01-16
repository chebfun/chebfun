% Test file for @chebfun/mtimes.m.

function pass = test_mtimes(pref)

% TODO: Test Chebfun2 outerproducts!

% Get preferences.
if (nargin < 1)
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Random numbers to use as arbitrary multiplicative constants.
alpha = -0.194758928283640 + 0.075474485412665i;

%% SCALAR-VALUED

% Check operation for empty inputs.
f = chebfun(@sin, [-1 1], pref);
pass(1) = isempty(f*[]);
pass(2) = isempty(f*chebfun());

% Turn on splitting, since we'll need it for the rest of the tests.
pref.splitting = 1;

% Test multiplication by scalars.
f1_op = @(x) sin(x).*abs(x - 0.1);
f1 = chebfun(f1_op, pref);
pass(3:4) = test_mult_function_by_scalar(f1, f1_op, alpha, x);

% Test use of mtimes() to compute products.
f = chebfun({@(x) sin(2*pi*x), @(x) cos(2*pi*x)}, [-1 0 1]);
g = chebfun({@(x) exp(2*pi*1i*x), @(x) cos(2*pi*x)}, [-1 0 1]);
tol = 10*max(f.vscale*f.epslevel, g.vscale*g.epslevel);
pass(5) = abs(g.'*f - (0.5 + 0.5i)) < tol;

%% ARRAY-VALUED

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

h = A*f2.';
h_exact = @(x) A*[sin(x).*abs(x - 0.1)  exp(x)].';
err = abs(feval(h, x) - h_exact(x));
pass(10) = max(err(:)) < 10*vscale(h)*epslevel(h);

%% QUASIMATRIX

% Test operation for array-valued chebfuns.
f2_op = @(x) [sin(x).*abs(x - 0.1)  exp(x)];
f2q = quasimatrix(f2_op, pref);
pass(11:12) = test_mult_function_by_scalar(f2q, f2_op, alpha, x);

f = chebfun({@(x) [sin(2*pi*x) sin(2*pi*x)], ...
    @(x) [cos(2*pi*x) cos(4*pi*x)]}, [-1 0 1]);
g = chebfun({@(x) [sin(4*pi*x) sin(4*pi*x)], ...
    @(x) [cos(2*pi*x) cos(4*pi*x)]}, [-1 0 1]);
fq = quasimatrix(f);
gq = quasimatrix(g);

err = gq.'*f - [0.5 0; 0 0.5];
tol = 10*max(vscale(f)*epslevel(f), vscale(g)*epslevel(g));
pass(13) = norm(err(:), inf) < tol;
err = g.'*fq - [0.5 0; 0 0.5];
tol = 10*max(vscale(f)*epslevel(f), vscale(g)*epslevel(g));
pass(14) = norm(err(:), inf) < tol;
err = gq.'*fq - [0.5 0; 0 0.5];
tol = 10*max(vscale(f)*epslevel(f), vscale(g)*epslevel(g));
pass(15) = norm(err(:), inf) < tol;

A = randn(2, 2);
g = f2q*A;
g_exact = @(x) f2_op(x)*A;
err = abs(feval(g, x) - g_exact(x));
pass(16) = max(err(:)) < 10*vscale(g)*epslevel(g);

h = A*f2q.';
h_exact = @(x) A*[sin(x).*abs(x - 0.1)  exp(x)].';
err = abs(feval(h, x) - h_exact(x));
pass(17) = max(err(:)) < 10*vscale(h)*epslevel(h);

%% Test error conditions.
try
    h = f*'X';
    pass(18) = false;
catch ME
    pass(18) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:mtimes:unknown');
end

try
    h = f*f1;
    pass(19) = false;
catch ME
    pass(19) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:mtimes:dims');
end

try
    h = f*g;
    pass(20) = false;
catch ME
    pass(20) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:mtimes:dims');
end

%% Test on SINGFUN - multiplication by scalar:

f = chebfun(@(x) sin(20*x)./((x+1).^0.5), 'exps', [-0.5 0], 'splitting', 'on');
h = 3*f;
h_op = @(x) 3*sin(20*x)./((x+1).^0.5);
h_vals = feval(h, x);
h_exact = h_op(x);
err = h_vals - h_exact;
pass(21) = norm(err, inf) < 2e2*vscale(h)*epslevel(h);

%% Test on SINGFUN - multiplication of a column CHEBFUN and a row CHEBFUN:CHEBFUN:

f = chebfun(@(x) sin(20*x)./((x+1).^0.5), 'exps', [-0.5 0], 'splitting', 'on');
f = f.';
g = chebfun(@(x) cos(30*x), 'splitting', 'on');
h = f*g;
h_exact = 0.13033807496531659;
err = h - h_exact;
pass(22) = abs(err) < 1e1*h_exact*max(epslevel(f), epslevel(g));

%% Tests for function defined on unbounded domain:

% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];
domCheck = [-1e6 -3*pi];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% Array-valued function:
op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
f = chebfun(op, dom);
A = randn(3, 3);
g = f*A;
gVals = feval(g, x);

op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x]*A;
gExact = op(x);
err = gVals - gExact;
pass(23) = norm(err, inf) < 1e1*max(get(g,'epslevel').*get(g,'vscale'));

% Some tests which make sure MTIMES behaves correctly regarding the dimensions
% of its inputs.
f = chebfun(@sin);
A = chebfun(@(x) [sin(x) cos(x) exp(x)]);
Q = cheb2quasi(A);

pass(24) = isequal(size(f'*f), [1 1]);
pass(25) = isequal(size(f*f'), [Inf Inf]);
pass(26) = isequal(size(f*[1 2]), [Inf 2]);
pass(27) = isequal(size([1 ; 2]*f'), [2 Inf]);
pass(28) = isequal(size(A'*A), [3 3]);
pass(29) = isequal(size(A*A'), [Inf Inf]);
pass(30) = isequal(size(f'*A), [1 3]);
pass(31) = isequal(size(A'*f), [3 1]);
pass(32) = isequal(size(A*ones(3, 2)), [Inf 2]);
pass(33) = isequal(size(A*ones(3, 3)), [Inf 3]);
pass(34) = isequal(size(A*ones(3, 4)), [Inf 4]);
pass(35) = isequal(size(ones(2, 3)*A'), [2 Inf]);
pass(36) = isequal(size(ones(3, 3)*A'), [3 Inf]);
pass(37) = isequal(size(ones(4, 3)*A'), [4 Inf]);
pass(38) = isequal(size(Q'*Q), [3 3]);
pass(39) = isequal(size(Q*Q'), [Inf Inf]);
pass(40) = isequal(size(f'*Q), [1 3]);
pass(41) = isequal(size(Q'*f), [3 1]);
pass(42) = isequal(size(Q*ones(3, 2)), [Inf 2]);
pass(43) = isequal(size(Q*ones(3, 3)), [Inf 3]);
pass(44) = isequal(size(Q*ones(3, 4)), [Inf 4]);
pass(45) = isequal(size(ones(2, 3)*Q'), [2 Inf]);
pass(46) = isequal(size(ones(3, 3)*Q'), [3 Inf]);
pass(47) = isequal(size(ones(4, 3)*Q'), [4 Inf]);
pass(48) = causesDimensionError(@() ones(3, 2)*A);
pass(49) = causesDimensionError(@() ones(3, 3)*A);
pass(50) = causesDimensionError(@() ones(3, 4)*A);
pass(51) = causesDimensionError(@() A'*ones(2, 3));
pass(52) = causesDimensionError(@() A'*ones(3, 3));
pass(53) = causesDimensionError(@() A'*ones(4, 3));
pass(54) = causesDimensionError(@() ones(3, 2)*Q);
pass(55) = causesDimensionError(@() ones(3, 3)*Q);
pass(56) = causesDimensionError(@() ones(3, 4)*Q);
pass(57) = causesDimensionError(@() Q'*ones(2, 3));
pass(58) = causesDimensionError(@() Q'*ones(3, 3));
pass(59) = causesDimensionError(@() Q'*ones(4, 3));
pass(60) = causesDimensionError(@() A*A);
pass(61) = causesDimensionError(@() A'*A');
pass(62) = causesDimensionError(@() Q*Q);
pass(63) = causesDimensionError(@() Q'*Q');

end

% Test the multiplication of a chebfun F, specified by F_OP, by a scalar ALPHA
% using a grid of points X in [-1  1] for testing samples.
function result = test_mult_function_by_scalar(f, f_op, alpha, x)
    g1 = f * alpha;
    g2 = alpha * f;
    result(1) = isequal(g1, g2);
    g_exact = @(x) f_op(x) * alpha;
    err = feval(g1, x) - g_exact(x);
    result(2) = norm(err(:), inf) < 10*max(vscale(g1).*epslevel(g1));
end

function result = causesDimensionError(op)
    try
        op();
        result = false;
    catch ME
        result = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:mtimes:dims');
    end
end

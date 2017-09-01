% Test file for @chebfun/times.m.

function pass = test_times(pref)

% Get preferences.
if ( nargin < 1 )
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
pass(1) = isempty(f.*[]);
pass(2) = isempty(f.*chebfun());

% Turn on splitting, since we'll need it for the rest of the tests.
pref.splitting = 1;

% Test multiplication by scalars.
f1_op = @(x) sin(x).*abs(x - 0.1);
f1 = chebfun(f1_op, pref);
pass(3:4) = test_mult_function_by_scalar(f1, f1_op, alpha, x);

% Test multiplication of two chebfun objects.
g1_op = @(x) cos(x).*sign(x + 0.2);
g1 = chebfun(g1_op, pref);
pass(5:6) = test_mult_function_by_function(f1, f1_op, g1, g1_op, x);

%% ARRAY-VALUED

% Test operation for array-valued chebfuns.
f2_op = @(x) [sin(x).*abs(x - 0.1)  exp(x)];
f2 = chebfun(f2_op, pref);
pass(7:8) = test_mult_function_by_scalar(f2, f2_op, alpha, x);

g2_op = @(x) [cos(x).*sign(x + 0.2) tan(x)];
g2 = chebfun(g2_op, pref);
pass(9:10) = test_mult_function_by_function(f2, f2_op, g2, g2_op, x);

% Test operation for transposed chebfuns.
pass(11:12) = test_mult_function_by_scalar(f1.', @(x) f1_op(x).', alpha, x);
pass(13:14) = test_mult_function_by_function(f1.', @(x) f1_op(x).', ...
    g1.', @(x) g1_op(x).', x);

%% QUASIMATRICES
f2q = quasimatrix(f2_op, pref);
pass(15:16) = test_mult_function_by_scalar(f2q, f2_op, alpha, x);

% Quasi * array-cheb
pass(17:18) = test_mult_function_by_function(f2q, f2_op, g2, g2_op, x);

% Quasi * quasi
g2q = quasimatrix(g2_op, pref);
pass(19:20) = test_mult_function_by_function(f2q, f2_op, g2q, g2_op, x);

%% Check error conditions.

try
    h = f1.*'X';
    pass(21) = false;
catch ME
    pass(21) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:times:unknown');
end

try
    h = f1.*f1.';
    pass(22) = false;
catch ME
    pass(22) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:times:matdim');
end

try
    h = f1.*g2q.';
    pass(23) = false;
catch ME
    pass(23) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:times:matdim');
end

%% Test on singular function:

dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

%% Case of a scalar and a function:
c = 3;
pow = -0.5;
op = @(x) (x - dom(2)).^pow.*sin(x);
op_exact = @(x) c*(x - dom(2)).^pow.*sin(x);
f = chebfun(op, dom, 'exps', [0 pow], 'splitting', 'on');
g = c.*f;
g_exact = chebfun(op_exact, dom, 'exps', [0 pow], 'splitting', 'on');

err = norm(feval(g, x) - feval(g_exact, x), inf);
pass(24) = ( err < 50*eps*norm(feval(g_exact, x), inf) );

%% Case of two functions: piecewise smooth chebfun - splitting on.
pow1 = -0.3;
pow2 = -0.5;
op1 = @(x) (x - dom(2)).^pow1.*sin(100*x);
op2 = @(x) (x - dom(2)).^pow2.*cos(300*x);
op_exact = @(x) (x - dom(2)).^(pow1+pow2).*sin(100*x).*cos(300*x);
f = chebfun(op1, dom, 'exps', [0 pow1], 'splitting', 'on');
g = chebfun(op2, dom, 'exps', [0 pow2] , 'splitting', 'on');
h = f.*g;
h_exact = chebfun(op_exact, dom, 'exps', [0 pow1+pow2], 'splitting', 'on');

err = norm(feval(h, x) - feval(h_exact, x), inf);
pass(25) = ( err < 1e4*eps*...
    norm(feval(h_exact, x), inf) );


%% Tests for function defined on unbounded domain:

% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];
domCheck = [-1e2 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

opf = @(x) x.^2.*exp(-x.^2);
opg = @(x) (1-exp(-x.^2))./x;
oph = @(x) x.*exp(-x.^2).*(1-exp(-x.^2));
f = chebfun(opf, dom);
g = chebfun(opg, dom);
h = f.*g;
hVals = feval(h, x);
hExact = oph(x);
err = hVals - hExact;
pass(26) = norm(err, inf) < 10*eps*get(f,'vscale');

%% Test multiplication between a CHEBFUN and a TRIGFUN.

dom = [0 pi 2*pi];

% 1. One column case.
f = chebfun(@(x) x.^2, dom, pref);
g = chebfun(@(x) cos(x), [dom(1) dom(end)], 'periodic');
h1 = f.*g;
% We want the result to use the same tech as the one used by f.
pass(27) = strcmpi(func2str(get(h1.funs{1}.onefun, 'tech')), ...
                   func2str(get(f.funs{1}.onefun, 'tech')));
h2 = chebfun(@(x) (x.^2).*cos(x), dom, pref);
pass(28) = norm(h1-h2, inf) < 1e1*eps*get(h2,'vscale');

% 2. Quasimatrix case.
f = chebfun(@(x) [cos(x), sin(x)], [dom(1) dom(end)], 'periodic');
g = chebfun(@(x) [x, x.^3], dom, pref);
h1 = f.*g;
% We want the result to use the same tech as the one used by g.
pass(29) = strcmpi(func2str(get(h1(:,1).funs{1}.onefun, 'tech')), ...
                   func2str(get(g(:,1).funs{1}.onefun, 'tech')));
pass(30) = strcmpi(func2str(get(h1(:,2).funs{1}.onefun, 'tech')), ...
                   func2str(get(g(:,2).funs{1}.onefun, 'tech')));
h2 = chebfun(@(x) [x.*cos(x), x.^3.*sin(x)], dom, pref);
pass(31) = norm(h1-h2, inf) < 1e2*eps*get(h2,'vscale');

end

%% The tests

% Test the multiplication of a chebfun F, specified by F_OP, by a scalar ALPHA
% using a grid of points X in [-1  1] for testing samples.
function result = test_mult_function_by_scalar(f, f_op, alpha, x)
    g1 = f .* alpha;
    g2 = alpha .* f;
    result(1) = isequal(g1, g2);
    g_exact = @(x) f_op(x) .* alpha;
    err = feval(g1, x) - g_exact(x);
    result(2) = norm(err(:), inf) < 1e2*max(vscale(g1)*eps);
        
end

% Test the multiplication of two chebfuns F and G, specified by F_OP and G_OP,
% using a grid of points X in [-1  1] for testing samples.
function result = test_mult_function_by_function(f, f_op, g, g_op, x)
    h1 = f .* g;
    h2 = g .* f;
    result(1) = norm(h1 - h2) < 10*max(vscale(h1)*eps);
    h_exact = @(x) f_op(x) .* g_op(x);
    err = feval(h1, x) - h_exact(x);
    result(2) = norm(err(:), inf) < 1e2*max(vscale(h1)*eps);
        
end

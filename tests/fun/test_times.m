% Test file for fun/times.m

function pass = test_times(pref)

% Get preferences.
if (nargin < 1)
    pref = chebpref();
end

% Set a domain
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(1000, 1) + dom(1);

% Random numbers to use as arbitrary multiplicative constants.
alpha = -0.194758928283640 + 0.075474485412665i;
beta = -0.526634844879922 - 0.685484380523668i;

%%
% Check operation in the face of empty arguments.

f = bndfun();
g = bndfun(@(x) x, dom, [], [], pref);
pass(1) = (isempty(f .* f) && isempty(f .* g) && isempty(g .* f));

%%
% Check multiplication by scalars.

f_op = @(x) sin(x);
f = bndfun(f_op, dom, [], [], pref);
pass(2:3) = test_mult_function_by_scalar(f, f_op, alpha, x);

f_op = @(x) [sin(x) cos(x)];
f = bndfun(f_op, dom, [], [], pref);
pass(4:5) = test_mult_function_by_scalar(f, f_op, alpha, x);

%%
% Check multiplication by constant functions.

f_op = @(x) sin(x);
f = bndfun(f_op, dom, [], [], pref);
g_op = @(x) alpha*ones(size(x));
g = bndfun(g_op, dom, [], [], pref);
pass(6) = test_mult_function_by_function(f, f_op, g, g_op, x, false);

% This should fail with a dimension mismatch error from bndfun.mtimes().
f_op = @(x) [sin(x) cos(x)];
f = bndfun(f_op, dom, [], [], pref);
g_op = @(x) repmat([alpha, beta], size(x, 1), 1);
g = bndfun(g_op, dom, [], [], pref);
pass(7) = test_mult_function_by_function(f, f_op, g, g_op, x, false);

%%
% Spot-check multiplication of two bndfun objects for a few test
% functions.

f_op = @(x) ones(size(x));
f = bndfun(f_op, dom, [], [], pref);
pass(8) = test_mult_function_by_function(f, f_op, f, f_op, x, false);

f_op = @(x) exp(x) - 1;
f = bndfun(f_op, dom, [], [], pref);

g_op = @(x) 1./(1 + x.^2);
g = bndfun(g_op, dom, [], [], pref);
pass(9) = test_mult_function_by_function(f, f_op, g, g_op, x, false);

g_op = @(x) cos(1e4*x);
g = bndfun(g_op, dom, [], [], pref);
pass(10) = test_mult_function_by_function(f, f_op, g, g_op, x, false);

g_op = @(t) sinh(t*exp(2*pi*1i/6));
g = bndfun(g_op, dom, [], [], pref);
pass(11) = test_mult_function_by_function(f, f_op, g, g_op, x, false);

%%
% Check operation for array-valued bndfun objects.

f = bndfun(@(x) [sin(x) cos(x) exp(x)], dom, [], [], pref);
g = bndfun(@(x) tanh(x), dom, [], [], pref);
h1 = f .* g;
h2 = g .* f;
pass(12) = isequal(h1, h2);
h_exact = @(x) [tanh(x).*sin(x) tanh(x).*cos(x) tanh(x).*exp(x)];
err = feval(h1, x) - h_exact(x);
pass(13) = max(abs(err(:))) < 10*max(get(h1, 'vscale').*get(h1, 'epslevel'));

g = bndfun(@(x) [sinh(x) cosh(x) tanh(x)], dom, [], [], pref);
h = f .* g;
h_exact = @(x) [sinh(x).*sin(x) cosh(x).*cos(x) tanh(x).*exp(x)];
err = feval(h, x) - h_exact(x);
pass(14) = max(abs(err(:))) < 10*max(get(h, 'vscale').*get(h, 'epslevel'));

% This should fail with a dimension mismatch error.
try
    g = bndfun(@(x) [sinh(x) cosh(x)], dom, [], [], pref);
    disp(f .* g);
    pass(15) = false;
catch ME
    pass(15) = strcmp(ME.identifier, 'CHEBFUN:CHEBTECH:times:dim2');
end

%%
% Check specially handled cases, including some in which an adjustment for
% positivity is performed.

f_op = @(t) sinh(t*exp(2*pi*1i/6));
f = bndfun(f_op, dom, [], [], pref);
pass(16) = test_mult_function_by_function(f, f_op, f, f_op, x, false);

g_op = @(t) conj(sinh(t*exp(2*pi*1i/6)));
g = conj(f);
pass(17:18) = test_mult_function_by_function(f, f_op, g, g_op, x, true);

f_op = @(x) exp(1i*x) - 1;
f = bndfun(f_op, dom, [], [], pref);
pass(19:20) = test_mult_function_by_function(f, f_op, f, f_op, x, false);

%%
% Check that multiplication and direct construction give similar results.

tol = 100*eps;
g_op = @(x) 1./(1 + x.^2);
g = bndfun(g_op, dom, [], [], pref);
h1 = f .* g;
h2 = bndfun(@(x) f_op(x) .* g_op(x), dom, [], [], pref);
pass(21) = normest(h1 - h2) < tol;

%%
% Check that multiplying a BNDFUN by an unhappy BNDFUN gives an unhappy
% result.

f = bndfun(@(x) cos(x+1), dom);    % Happy
g = bndfun(@(x) sqrt(x+1), dom);   % Unhappy
h = f.*g;  % Multiply unhappy by happy.
pass(22) = (~get(g, 'ishappy')) && (~get(h, 'ishappy')); %#ok<*BDSCI,*BDLGI>
h = g.*f;  % Multiply happy by unhappy.
pass(23) = (~get(g, 'ishappy')) && (~get(h, 'ishappy'));

end

% Test the multiplication of a BNDFUN F, specified by F_OP, by a scalar ALPHA
% using a grid of points X in [a  b] for testing samples.
function result = test_mult_function_by_scalar(f, f_op, alpha, x)
g1 = f .* alpha;
g2 = alpha .* f;
result(1) = isequal(g1, g2);
g_exact = @(x) f_op(x) .* alpha;
tol = 10*max(get(g1, 'vscale').*get(g1, 'epslevel'));
result(2) = norm(feval(g1, x) - g_exact(x), inf) < tol;
end

% Test the multiplication of two BNDFUN objects F and G, specified by F_OP and
% G_OP, using a grid of points X in [a  b] for testing samples.  If CHECKPOS is
% TRUE, an additional check is performed to ensure that the values of the result
% are all nonnegative; otherwise, this check is skipped.
function result = test_mult_function_by_function(f, f_op, g, g_op, x, checkpos)
h = f .* g;
h_exact = @(x) f_op(x) .* g_op(x);
tol = 10*max(get(h, 'vscale').*get(h, 'epslevel'));
result(1) = all(max(abs(feval(h, x) - h_exact(x))) < 5*tol);
if ( checkpos )
    result(2) = all(feval(h, x) >= 0);
end

end

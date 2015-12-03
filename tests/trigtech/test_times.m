% Test file for trigtech/times.m

function pass = test_times(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Random numbers to use as arbitrary multiplicative constants.
alpha = -0.194758928283640 + 0.075474485412665i;
beta = -0.526634844879922 - 0.685484380523668i;

%%
% Check operation in the face of empty arguments.

f = testclass.make();
g = testclass.make(@(x) sin(pi*x), [], pref);
pass(1) = (isempty(f .* f) && isempty(f .* g) && isempty(g .* f));

%%
% Check multiplication by scalars.

f_op = @(x) sin(cos(pi*x));
f = testclass.make(f_op, [], pref);
pass(2:3) = test_mult_function_by_scalar(f, f_op, alpha, x);

f_op = @(x) exp([sin(pi*x) -cos(pi*x)]);
f = testclass.make(f_op, [], pref);
pass(4:5) = test_mult_function_by_scalar(f, f_op, alpha, x);

%%
% Check multiplication by constant functions.

f_op = @(x) 3./(4-cos(pi*x));
f = testclass.make(f_op, [], pref);
g_op = @(x) alpha*ones(size(x));
g = testclass.make(g_op, [], pref);
pass(6) = test_mult_function_by_function(f, f_op, g, g_op, x, false);

% This should fail with a dimension mismatch error from trigtech.mtimes().
f_op = @(x) [sin(pi*x) cos(pi*x)];
f = testclass.make(f_op, [], pref);
g_op = @(x) repmat([alpha, beta], size(x, 1), 1);
g = testclass.make(g_op, [], pref);
pass(7) = test_mult_function_by_function(f, f_op, g, g_op, x, false);

%%
% Spot-check multiplication of two trigtech objects for a few test
% functions.

f_op = @(x) ones(size(x));
f = testclass.make(f_op, [], pref);
pass(8) = test_mult_function_by_function(f, f_op, f, f_op, x, false);

f_op = @(x) exp(cos(pi*x)) - 1;
f = testclass.make(f_op, [], pref);

g_op = @(x) 3./(4-cos(pi*x));
g = testclass.make(g_op, [], pref);
pass(9) = test_mult_function_by_function(f, f_op, g, g_op, x, false);

g_op = @(x) cos(1e4*pi*x);
g = testclass.make(g_op, [], pref);
pass(10) = test_mult_function_by_function(f, f_op, g, g_op, x, false);

g_op = @(x) exp(1i*1e2*pi*x);
g = testclass.make(g_op, [], pref);
pass(11) = test_mult_function_by_function(f, f_op, g, g_op, x, false);

%%
% Check operation for array-valued trigtech objects.
f_op = @(x) [sin(pi*x) cos(30*pi*x) 3./(4-cos(pi*x))];
f = testclass.make(f_op, [], pref);
g_op = @(x) tanh(sin(pi*x)+cos(pi*x));
g = testclass.make(g_op, [], pref);
h1 = f .* g;
h2 = g .* f;
pass(12) = (norm(h1.coeffs-h2.coeffs) < 10*eps);
h_exact = @(x) bsxfun(@times,g_op(x),f_op(x));
err = feval(h1, x) - h_exact(x);
pass(13) = max(abs(err(:))) < 100*eps;
        
g_op = @(x) [tanh(sin(pi*x)+cos(pi*x)) sin(pi*x) exp(sin(pi*x))];
g = testclass.make(g_op, [], pref);
h = f .* g;
h_exact = @(x) g_op(x).*f_op(x);
err = feval(h, x) - h_exact(x);
pass(14) = max(abs(err(:))) < 100*eps;

% This should fail with a dimension mismatch error.
try
    g = testclass.make(@(x) [sin(pi*x) cos(pi*x)], [], pref);
    disp(f .* g);
    pass(15) = false;
catch ME
    pass(15) = strcmp(ME.identifier, 'CHEBFUN:TRIGTECH:times:dim2');
end

%%
% Check specially handled cases, including some in which an adjustment for
% positivity is performed.

f_op = @(x) exp(cos(pi*x)) + exp(1i*2*pi*x);
f = testclass.make(f_op, [], pref);
pass(16) = test_mult_function_by_function(f, f_op, f, f_op, x, false);

g_op = @(t) conj(exp(cos(pi*x)) + exp(1i*2*pi*x));
g = conj(f);
pass(17:18) = test_mult_function_by_function(f, f_op, g, g_op, x, true);

f_op = @(x) 1+cos(pi*x);
f = testclass.make(f_op, [], pref);
pass(19:20) = test_mult_function_by_function(f, f_op, f, f_op, x, true);

%%
% Check that multiplication and direct construction give similar results.

tol = 50*eps;
g_op = @(x) 3./(4 - cos(2*pi*x));
g = testclass.make(g_op, [], pref);
h1 = f .* g;
h2 = testclass.make(@(x) f_op(x) .* g_op(x), [], pref);
h2 = prolong(h2, length(h1));
pass(21) = norm(h1.coeffs - h2.coeffs, inf) < tol;

%%
% Check that multiplying a TRIGTECH by an unhappy TRIGTECH gives an unhappy
% result.

warning off; % Suppress expected warnings about unhappy operations.
f = testclass.make(@(x) exp(cos(pi*x)));    % Happy
g = testclass.make(@(x) x);   % Unhappy
h = f.*g;  % Multiply unhappy by happy.
pass(22) = (~g.ishappy) && (~h.ishappy); %#ok<*BDSCI,*BDLGI>
h = g.*f;  % Multiply happy by unhappy.
pass(23) = (~g.ishappy) && (~h.ishappy);
warning on; % Re-enable warnings.
end

% Test the multiplication of a TRIGTECH F, specified by F_OP, by a scalar ALPHA
% using a grid of points X in [-1  1] for testing samples.
function result = test_mult_function_by_scalar(f, f_op, alpha, x)
    g1 = f .* alpha;
    g2 = alpha .* f;
    result(1) = isequal(g1, g2);
    g_exact = @(x) f_op(x) .* alpha;
    result(2) = norm(feval(g1, x) - g_exact(x), inf) < ...
        200*max(vscale(g1)*eps);
end

% Test the multiplication of two TRIGTECH objects F and G, specified by F_OP and
% G_OP, using a grid of points X in [-1  1] for testing samples.  If CHECKPOS is
% TRUE, an additional check is performed to ensure that the values of the result
% are all nonnegative; otherwise, this check is skipped.
function result = test_mult_function_by_function(f, f_op, g, g_op, x, checkpos)
    h = f .* g;
    h_exact = @(x) f_op(x) .* g_op(x);
    result(1) = norm(feval(h, x) - h_exact(x), inf) < ...
        1e5*max(vscale(h)*eps);
        
    if ( checkpos )
        values = h.coeffs2vals(h.coeffs); 
        result(2) = all(values >= 0);
    end
end

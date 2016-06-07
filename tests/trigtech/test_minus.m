% Test file for trigtech/minus.m

function pass = test_minus(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% A random number to use as an arbitrary additive constant.
alpha = -0.194758928283640 + 0.075474485412665i;

%%
% Check operation in the face of empty arguments.

f = testclass.make();
g = testclass.make(@(x) cos(pi*x), [], pref);
pass(1) = (isempty(f - f) && isempty(f - g) && isempty(g - f));

%%
% Check subtraction with scalars.

f_op = @(x) exp(sin(pi*x));
f = testclass.make(f_op, [], pref);
pass(2:3) = test_sub_function_and_scalar(f, f_op, alpha, x);

%%
% Check subtraction of two TRIGTECH objects.

f_op = @(x) zeros(size(x));
f = testclass.make(f_op, [], pref);
pass(4:5) = test_sub_function_and_function(f, f_op, f, f_op, x);

f_op = @(x) exp(cos(pi*x)) - 1;
f = testclass.make(f_op, [], pref);

g_op = @(x) sin(100*pi*x);
g = testclass.make(g_op, [], pref);
pass(6:7) = test_sub_function_and_function(f, f_op, g, g_op, x);

g_op = @(x) sin(cos(10*pi*x));
g = testclass.make(g_op, [], pref);
pass(8:9) = test_sub_function_and_function(f, f_op, g, g_op, x);
    
%%
% Check operation for array-valued TRIGTECH objects.

f_op = @(x) [zeros(size(x)) zeros(size(x)) zeros(size(x))];
f = testclass.make(f_op, [], pref);
pass(10:11) = test_sub_function_and_function(f, f_op, f, f_op, x);

f_op = @(x) [sin(10*pi*x) sin(cos(pi*x)) exp(cos(pi*x))];
f = testclass.make(f_op, [], pref);
pass(12:13) = test_sub_function_and_scalar(f, f_op, alpha, x);

g_op = @(x) [sin(pi*x) exp(1i*pi*x).*exp(1i*pi*x) cos(pi*x)];
g = testclass.make(g_op, [], pref);
pass(14:15) = test_sub_function_and_function(f, f_op, g, g_op, x);

% This should fail with a dimension mismatch error.
g_op = @(x) sin(10*pi*x);
g = testclass.make(g_op, [], pref);
try
    h = f - g; %#ok<NASGU>
    pass(16) = false;
catch ME
    pass(16) = strcmp(ME.message, 'Matrix dimensions must agree.');
end

%%
% Check that direct construction and MINUS give comparable results.

tol = 10*eps;
f = testclass.make(@(x) sin(pi*cos(3*pi*x)), [], pref);
g = testclass.make(@(x) cos(pi*sin(10*pi*x)) - 1, [], pref);
h1 = f - g;
h2 = testclass.make(@(x) sin(pi*cos(3*pi*x)) - (cos(pi*sin(10*pi*x)) - 1), [], pref);
pass(17) = norm(h1.coeffs - h2.coeffs, inf) < tol;

%%
% Check that subtracting a TRIGTECH and an unhappy TRIGTECH gives an
% unhappy result.

f = testclass.make(@(x) cos(pi*x));    % Happy
g = testclass.make(@(x) cos(x));   % Unhappy
h = f - g;  % Subtract unhappy from happy.
pass(18) = (~g.ishappy) && (~h.ishappy);
h = g - f;  % Subtract happy from unhappy.
pass(19) = (~g.ishappy) && (~h.ishappy);

%%
% Test subtraction of array-valued scalar to array-valued TRIGTECH.

f = testclass.make(@(x) exp([sin(pi*x) cos(pi*x) -sin(pi*x).^2]));
g = f - [1 2 3];
g_exact = @(x) [exp(sin(pi*x))-1 exp(cos(pi*x))-2 exp(-sin(pi*x).^2)-3];
err = feval(g, x) - g_exact(x);
pass(20) = norm(err(:), inf) < 10*max(vscale(g)*eps);

%%
% Test scalar expansion in TRIGTECH argument.

f = testclass.make(@(x) sin(pi*x));
g = f - [1 2 3];
g_exact = @(x) [(-1 + sin(pi*x)) (-2 + sin(pi*x)) (-3 + sin(pi*x))];
err = feval(g, x) - g_exact(x);
pass(21) = isequal(size(g.coeffs, 2), 3) && norm(err(:), inf) < ...
    10*max(vscale(g)*eps);

end

% Test the subtraction of a TRIGTECH F, specified by F_OP, to and from a scalar
% ALPHA using a grid of points X in [-1  1] for testing samples.
function result = test_sub_function_and_scalar(f, f_op, alpha, x)
    g1 = f - alpha;
    g2 = alpha - f;
    result(1) = isequal(g1, -g2);
    g_exact = @(x) f_op(x) - alpha;
    result(2) = norm(feval(g1, x) - g_exact(x), inf) <= ...
        1000*max(vscale(g1)*eps);
end

% Test the subraction of two TRIGTECH objects F and G, specified by F_OP and
% G_OP, using a grid of points X in [-1  1] for testing samples.
function result = test_sub_function_and_function(f, f_op, g, g_op, x)
    h1 = f - g;
    h2 = g - f;
    result(1) = isequal(h1, -h2);
    h_exact = @(x) f_op(x) - g_op(x);
    norm(feval(h1, x) - h_exact(x), inf);
    result(2) = norm(feval(h1, x) - h_exact(x), inf) <= ...
        1e3*max(vscale(h1)*eps);       
end

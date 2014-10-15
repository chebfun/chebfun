% Test file for trigtech/rdivide.m

function pass = test_rdivide(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(10, 1) - 1;

% Random numbers to use as arbitrary constants.
alpha = -0.194758928283640 + 0.075474485412665i;
beta = -0.526634844879922 - 0.685484380523668i;

%%
% Check division by single scalars.

f_op = @(x) sin(10*pi*x);
f = testclass.make(f_op, [], pref);
pass(1) = test_div_function_by_scalar(f, f_op, alpha, x);

g = f ./ 0;
pass(2) = isnan(g);

f_op = @(x) [sin(10*pi*x) sin(20*pi*x)];
f = testclass.make(f_op, [], pref);
pass(3) = test_div_function_by_scalar(f, f_op, alpha, x);

g = f ./ 0;
pass(4) = isnan(g);

%%
% Check division by a row matrix of scalars in the case of an array-valued
% TRIGTECH object.

g = f ./ [alpha beta];
g_exact = @(x) [sin(10*pi*x)./alpha sin(20*pi*x)./beta];
pass(5) = norm(feval(g, x) - g_exact(x), inf) < 100*max(g.epslevel);

g = f ./ [alpha 0];
pass(6) = isnan(g) && ~any(isnan(g.coeffs(:, 1))) ...
    && all(isnan(g.coeffs(:, 2)));

%%
% Check division of a scalar by a TRIGTECH object.

f_op = @(x) exp(cos(pi*x));
f = testclass.make(@(x) exp(cos(pi*x)), [], pref);
pass(7) = test_div_scalar_by_function(alpha, f, f_op, x);

%%
% Check division of two TRIGTECH objects.

g_op = @(x) exp(cos(20*pi*x));
g = testclass.make(g_op, [], pref);

f_op = @(x) exp(cos(20*pi*x)) - 1;
f = testclass.make(f_op, [], pref);
pass(8) = test_div_function_by_function(f, f_op, g, g_op, x);

f_op = @(x) cos(1e3*pi*x);
f = testclass.make(f_op, [], pref);
pass(9) = test_div_function_by_function(f, f_op, g, g_op, x);

%%
% Check proper behavior under error conditions.

% Can't divide by a scalar matrix with multiple rows.
try
    f = testclass.make(@(x) sin(x), [], pref);
    disp(f ./ [1 ; 2]);
    pass(10) = false;
catch ME
    pass(10) = strcmp(ME.identifier, 'CHEBFUN:TRIGTECH:rdivide:size');
end

% Can't divide by a scalar row matrix if the column counts don't match.
try
    f = testclass.make(@(x) [sin(x) cos(x)], [], pref);
    disp(f ./ [1 2 3]);
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.identifier, 'CHEBFUN:TRIGTECH:rdivide:size');
end

%%
% Check that direct construction and RDIVIDE give comparable results.

tol = 100*eps;

xx = linspace(-1,1,10);
f = testclass.make(@(x) sin(10*pi*x), [], pref);
h1 = f ./ alpha;
h2 = testclass.make(@(x) sin(10*pi*x) ./ alpha, [], pref);
pass(12) = norm(feval(h1, xx) - feval(h2,xx), inf) < tol;

g = testclass.make(@(x) exp(cos(pi*x)), [], pref);
h1 = f ./ g;
h2 = testclass.make(@(x) sin(10*pi*x) ./ exp(cos(pi*x)), [], pref);
pass(13) = norm(feval(h1, xx) - feval(h2,xx), inf) < tol;
end

% Test the division of a TRIGTECH F, specified by F_OP, by a scalar ALPHA using
% a grid of points X in [-1  1] for testing samples.
function result = test_div_function_by_scalar(f, f_op, alpha, x)
    g = f ./ alpha;
    g_exact = @(x) f_op(x) ./ alpha;
    result = norm(feval(g, x) - g_exact(x), inf) < 10*max(g.vscale.*g.epslevel);
end

% Test the division of a scalar ALPHA by a TRIGTECH, specified by F_OP, using
% a grid of points X in [-1  1] for testing samples.
function result = test_div_scalar_by_function(alpha, f, f_op, x)
    g = alpha ./ f;
    g_exact = @(x) alpha ./ f_op(x);
    err = norm(feval(g, x) - g_exact(x), inf);
    result = err < 100*max(g.vscale.*g.epslevel);
end

% Test the division of two TRIGTECH objects F and G, specified by F_OP and
% G_OP, using a grid of points X in [-1  1] for testing samples.
function result = test_div_function_by_function(f, f_op, g, g_op, x)
    h = f ./ g;
    h_exact = @(x) f_op(x) ./ g_op(x);
    err = norm(feval(h, x) - h_exact(x), inf);
    result = err < 100*max(h.vscale.*h.epslevel);
end

% Test file for funcheb2/rdivide.

function pass = test_rdivide(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb2.pref();
end

% Generate a few random points to use as test values.
rand('seed', 6178); %#ok<RAND>
x = 2 * rand(100, 1) - 1;

% Random numbers to use as arbitrary constants.
alpha = -0.194758928283640 + 0.075474485412665i;
beta = -0.526634844879922 - 0.685484380523668i;

%%
% Check division by single scalars.

f_op = @(x) sin(x);
f = funcheb2(f_op, pref);
pass(1) = test_div_function_by_scalar(f, f_op, alpha, x);

g = f ./ 0;
pass(2) = isnan(g);

f_op = @(x) [sin(x) cos(x)];
f = funcheb2(f_op, pref);
pass(3) = test_div_function_by_scalar(f, f_op, alpha, x);

g = f ./ 0;
pass(4) = isnan(g);

%%
% Check division by a row matrix of scalars in the case of a vectorized
% funcheb2 object.

g = f ./ [alpha beta];
g_exact = @(x) [sin(x)./alpha cos(x)./beta];
pass(5) = norm(feval(g, x) - g_exact(x), 'inf') < 10*g.epslevel;

g = f ./ [alpha 0];
pass(6) = isnan(g) && ~any(isnan(g.values(:, 1))) ...
                   && all(isnan(g.values(:, 2)));

%%
% Check division of a scalar by a funcheb2 object.

f_op = @(x) exp(x);
f = funcheb2(@(x) exp(x), pref);
pass(7) = test_div_scalar_by_function(alpha, f, f_op, x);

%%
% Check division of two funcheb2 objects.

g_op = @(x) exp(x);
g = funcheb2(g_op, pref);

f_op = @(x) exp(x) - 1;
f = funcheb2(f_op, pref);
pass(8) = test_div_function_by_function(f, f_op, g, g_op, x);

f_op = @(x) 1./(1 + x.^2);
f = funcheb2(f_op, pref);
pass(9) = test_div_function_by_function(f, f_op, g, g_op, x);

f_op = @(x) cos(1e4*x);
f = funcheb2(f_op, pref);
pass(10) = test_div_function_by_function(f, f_op, g, g_op, x);

f_op = @(t) sinh(t*exp(2*pi*1i/6));
f = funcheb2(f_op, pref);
pass(11) = test_div_function_by_function(f, f_op, g, g_op, x);

%%
% Check proper behavior under error conditions.

% Can't divide by a scalar matrix with multiple rows.
try
    f = funcheb2(@(x) sin(x), pref);
    disp(f ./ [1 ; 2]);
    pass(12) = false;
catch ME
    pass(12) = strcmp(ME.identifier, 'CHEBUFN:FUNCHEB:rdivide:size');
end

% Can't divide by a scalar row matrix if the column counts don't match.
try
    f = funcheb2(@(x) [sin(x) cos(x)], pref);
    disp(f ./ [1 2 3]);
    pass(13) = false;
catch ME
    pass(13) = strcmp(ME.identifier, 'CHEBUFN:FUNCHEB:rdivide:size');
end

% Can't divide by a function which has roots inside [-1, 1].
try
    f = funcheb2(@(x) exp(x));
    g = funcheb2(@(x) sin(x));
    disp(f ./ g);
    pass(14) = false;
catch ME
    pass(14) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:rdivide:DivideByZeros');
end

%%
% Check that direct construction and RDIVIDE give comparable results.

tol = 100*eps;

f = funcheb2(@(x) sin(x), pref);
h1 = f ./ alpha;
h2 = funcheb2(@(x) sin(x) ./ alpha, pref);
pass(15) = norm(h1.values - h2.values, 'inf') < tol;

g = funcheb2(@(x) exp(x), pref);
h1 = f ./ g;
h2 = funcheb2(@(x) sin(x) ./ exp(x), pref);
pass(16) = norm(h1.values - h2.values, 'inf') < tol;

end

% Test the division of a FUNCHEB2 F, specified by F_OP, by a scalar ALPHA using
% a grid of points X in [-1  1] for testing samples.
function result = test_div_function_by_scalar(f, f_op, alpha, x)
    g = f ./ alpha;
    g_exact = @(x) f_op(x) ./ alpha;
    result = norm(feval(g, x) - g_exact(x), 'inf') < 10*g.epslevel;
end

% Test the division of a scalar ALPHA by a FUNCHEB2, specified by F_OP, using
% a grid of points X in [-1  1] for testing samples.
function result = test_div_scalar_by_function(alpha, f, f_op, x)
    g = alpha ./ f;
    g_exact = @(x) alpha ./ f_op(x);
    result = norm(feval(g, x) - g_exact(x), 'inf') < 10*g.epslevel;
end

% Test the division of two FUNCHEB2 objects F and G, specified by F_OP and
% G_OP, using a grid of points X in [-1  1] for testing samples.
function result = test_div_function_by_function(f, f_op, g, g_op, x)
    h = f ./ g;
    h_exact = @(x) f_op(x) ./ g_op(x);
    norm(feval(h, x) - h_exact(x), 'inf');
    result = norm(feval(h, x) - h_exact(x), 'inf') < 10*h.epslevel;
end

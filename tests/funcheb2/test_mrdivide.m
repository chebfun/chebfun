% Test file for funcheb2/mrdivide.

function pass = test_mrdivide(pref)

% Get preferences.
if (nargin < 1)
    pref = funcheb2.pref();
end

% Generate a few random points to use as test values.
rngstate = rng();
rng(6178);
x = 2 * rand(100, 1) - 1;

% Random number to use as a scalar constant.
alpha = -0.194758928283640 + 0.075474485412665i;

%%
% Check division of a funcheb2 object by a numeric array.

f_op = @(x) [sin(x) cos(x)];
f = funcheb2(f_op, pref);
pass(1) = isnan(f / 0);

g = f / alpha;
g_exact = @(x) f_op(x) ./ alpha;
pass(2) = norm(feval(g, x) - g_exact(x), 'inf') < 10*g.epslevel;

% A "least-squares" case where the solution is obvious.
I = eye(2);
g = f / I;
err = g*I - f;
pass(3) = max(max(abs(feval(err, x)))) < 10*g.epslevel;

% A less trivial least-squares case for which we still know the answer.
A = [1 1];
g = f / A;
g_exact = @(x) (sin(x) + cos(x))/2;
pass(4) = norm(feval(g, x) - g_exact(x), 'inf') < 10*g.epslevel;

%%
% Check division of a numeric array by a funcheb2 object.

f = funcheb2(@(x) sin(x));
g = alpha / f;
pass(5) = abs(innerProduct(f, g) - alpha) < 10*g.epslevel;

f = funcheb2(@(x) [sin(2*pi*x) cos(2*pi*x)]);
g = [1 1]/f;
g_exact = @(x) (sin(2*pi*x) + cos(2*pi*x));
pass(6) = norm(feval(g, x) - g_exact(x), 'inf') < 10*g.epslevel;

%%
% Check error conditions.

% Catch dimension mismatch errors.
try
    g = f / [1 2 3];
    pass(7) = false;
catch ME
    pass(7) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:mrdivide:size');
end

% Can't do f/g if both f and g are funcheb2 objects.
try
    f = funcheb2(@(x) sin(x));
    g = funcheb2(@(x) cos(x));
    h = f / g;
    pass(8) = false;
catch ME
    pass(8) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:mrdivide:funfun');
end

% Can't call mldivide on a funcheb2 and a non-funcheb2 or non-double object.
try
    f = funcheb2(@(x) sin(x));
    g = f / true;
    pass(9) = false;
catch ME
    pass(9) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:mrdivide:derp');
end

%%
% Restore the RNG state.

rng(rngstate);

end


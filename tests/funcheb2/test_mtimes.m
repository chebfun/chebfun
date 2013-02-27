% Test file for funcheb2/mtimes.

function pass = test_mtimes(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb2.pref();
end

% Generate a few random points to use as test values.
rngstate = rng();
rng(6178);
x = 2 * rand(100, 1) - 1;

% A random number to use as an arbitrary scalar multiplier.
alpha = randn() + 1i*randn();

%%
% Check operation in the face of empty arguments.

f = funcheb2(@(x) sin(x), pref);
g = funcheb2();
pass(1) = isempty(f*[]) && isempty([]*f) && isempty(2*g) && isempty(g*2);

%%
% Check operation for scalar funcheb2 objects.

f = funcheb2(@(x) sin(x), pref);
g1 = alpha*f;
g2 = f*alpha;
pass(2) = isequal(g1, g2);
g_exact = @(x) alpha*sin(x);
pass(3) = norm(feval(g1, x) - g_exact(x), 'inf') < 10*g1.epslevel;

g = 0*f;
pass(4) = all(g.values == 0) && all(g.coeffs == 0);

%%
% Check operation for vectorized funcheb2 objects.

f = funcheb2(@(x) [sin(x) cos(x) exp(x)], pref);
g1 = alpha*f;
g2 = f*alpha;
pass(5) = isequal(g1, g2);
g_exact = @(x) alpha*[sin(x) cos(x) exp(x)];
err = abs(feval(g1, x) - g_exact(x));
pass(6) = max(err(:)) < 10*g1.epslevel;

g = 0*f;
pass(7) = all(g.values == 0) && all(g.coeffs == 0);

A = randn(3, 3);
g = f*A;
g_exact = @(x) [sin(x) cos(x) exp(x)]*A;
err = abs(feval(g, x) - g_exact(x));
pass(8) = max(err(:)) < 10*g.epslevel;

%%
% Verify error handling and corner cases.

% Multiply non-scalar double and funcheb2.
try
    f = funcheb2(@(x) exp(x));
    disp([1 2 3]*f)
catch ME
    pass(9) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:mtimes:size') ...
        && strcmp(ME.message, 'Inner matrix dimensions must agree.');
end

% Multiply funcheb2 and non-sclar double with mismatching dimensions.
try
    f = funcheb2(@(x) [sin(x) cos(x)]);
    disp(f*[1 ; 2 ; 3]);
    pass(10) = false;
catch ME
    pass(10) = strcmp(ME.identifier, 'CHEBFUN:FUNCHEB:mtimes:size2') ...
        && strcmp(ME.message, 'Inner matrix dimensions must agree.');
end

% Using * for multiplication of two funcheb2 objects.
try
    g = funcheb2(@(x) x);
    disp(f*g);
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.message, 'Use .* to multiply FUNCHEB objects.');
end

% Using * to multiply a funcheb2 and something else.
try
    disp(f*uint8(128));
    pass(12) = false;
catch ME
    pass(12) = strcmp(ME.message, ...
        'mtimes does not know how to multiply a FUNCHEB and a uint8.');
end

%%
% Restore the RNG state.

rng(rngstate);

end


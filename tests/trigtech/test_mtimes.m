% Test file for trigtechh/mtimes.m

function pass = test_mtimes(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% A random number to use as an arbitrary scalar multiplier.
alpha = randn() + 1i*randn();

%%
% Check operation in the face of empty arguments.

f = testclass.make(@(x) sin(10*pi*x), [], pref);
g = testclass.make();
pass(1) = isempty(f*[]) && isempty([]*f) && isempty(2*g) && isempty(g*2);

%%
% Check operation for scalar TRIGTECH objects.

f = testclass.make(@(x) sin(10*pi*x), [], pref);
g1 = alpha*f;
g2 = f*alpha;
pass(2) = isequal(g1, g2);
g_exact = @(x) alpha*sin(10*pi*x);
pass(3) = norm(feval(g1, x) - g_exact(x), inf) < ...
    10*g1.vscale.*g1.epslevel;

g = 0*f;
pass(4) = all(g.coeffs == 0);

%%
% Check operation for array-valued TRIGTECH objects.

f = testclass.make(@(x) [sin(10*pi*x) cos(20*pi*x) cos(sin(pi*x))], [], pref);
g1 = alpha*f;
g2 = f*alpha;
pass(5) = isequal(g1, g2);
g_exact = @(x) alpha*[sin(10*pi*x) cos(20*pi*x) cos(sin(pi*x))];
err = abs(feval(g1, x) - g_exact(x));
pass(6) = max(err(:)) < 10*max(g1.vscale.*g1.epslevel);

g = 0*f;
pass(7) =  all(g.coeffs == 0);

A = randn(3, 3);
g = f*A;
g_exact = @(x) [sin(10*pi*x) cos(20*pi*x) cos(sin(pi*x))]*A;
err = abs(feval(g, x) - g_exact(x));
pass(8) = max(err(:)) < 10*max(g.vscale.*g.epslevel);

f = testclass.make(@(x) [exp(1i*11*pi*x) cos(20*pi*x) cos(sin(pi*x))], [], pref);
g = f*A;
g_exact = @(x) [exp(1i*11*pi*x) cos(20*pi*x) cos(sin(pi*x))]*A;
err = abs(feval(g, x) - g_exact(x));
pass(9) = max(err(:)) < 10*max(g.vscale.*g.epslevel);

f = testclass.make(@(x) [exp(1i*11*pi*x) cos(20*pi*x) cos(sin(pi*x))], [], pref);
A = randn(3, 3) + 1i*randn(3, 3);
g = f*A;
g_exact = @(x) [exp(1i*11*pi*x) cos(20*pi*x) cos(sin(pi*x))]*A;
err = abs(feval(g, x) - g_exact(x));
pass(10) = max(err(:)) < 10*max(g.vscale.*g.epslevel);

%%
% Verify error handling and corner cases.

% Multiply non-scalar double and TRIGTECH.
try
    f = testclass.make(@(x) exp(cos(pi*x)));
    disp([1 2 3]*f);
    pass(11) = false;
catch ME
    pass(11) = strcmp(ME.identifier, 'CHEBFUN:TRIGTECH:mtimes:size') ...
        && strcmp(ME.message, 'Inner matrix dimensions must agree.');
end

% Multiply TRIGTECH and non-scalar double with mismatching dimensions.
try
    f = testclass.make(@(x) [sin(10*pi*x) cos(20*pi*x)]);
    disp(f*[1 ; 2 ; 3]);
    pass(12) = false;
catch ME
    pass(12) = strcmp(ME.identifier, 'CHEBFUN:TRIGTECH:mtimes:size2') ...
        && strcmp(ME.message, 'Inner matrix dimensions must agree.');
end

% Using * for multiplication of two TRIGTECH objects.
try
    g = testclass.make(@(x) cos(20*pi*x));
    disp(f*g);
    pass(13) = false;
catch ME
    pass(13) = strcmp(ME.message, 'Use .* to multiply TRIGTECH objects.');
end

% Using * to multiply a TRIGTECH and something else.
try
    disp(f*uint8(128));
    pass(14) = false;
catch ME
    pass(14) = strcmp(ME.message, ...
        'mtimes does not know how to multiply a TRIGTECH and a uint8.');
end

end

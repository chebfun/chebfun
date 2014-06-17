% Test file for @chebfun/unwrap.m.

function pass = test_unwrap(pref)

if (nargin < 1)
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Test empty case.
pass(1) = isempty(unwrap(chebfun()));

% Check that chebfuns with only one fun aren't unwrapped.
f = chebfun(@(x) exp(x), [-1 1], pref);
pass(2) = isequal(f, unwrap(f));

pref.techPrefs.extrapolate = 1;
f = chebfun(@mysawtooth, [0 2*pi 4*pi 6*pi 8*pi], pref);
uf = unwrap(f);
pass(3) = norm(feval(uf, x) - (x - pi), inf) < 10*vscale(uf)*epslevel(uf);

% Check tol input.
g = 10*f;
ug = unwrap(g, 10*pi);
pass(4) = norm(feval(ug, x) - 10*(x - pi), inf) < 10*vscale(ug)*epslevel(ug);

% Check error conditions.
try
    f = chebfun(@(x) [sin(x) cos(x)], [-1 0 1]);
    g = unwrap(f);
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:unwrap:array');
end

end

function y = mysawtooth(x)
    y = mod(x, 2*pi) - pi;
end

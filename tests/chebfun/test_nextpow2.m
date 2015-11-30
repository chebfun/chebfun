% Test file for @chebfun/nextpow2.m.

function pass = test_nextpow2(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% NEXTPOW2 is not vectorized in versions of MATLAB prior to R2010a.  Vectorize
% it manually if we're running on older platforms.
if ( verLessThan('matlab', '7.10') )
    mynextpow2 = @nextpow2Vectorized;
else
    mynextpow2 = @nextpow2;
end

% Test for scalar-valued chebfuns.
f = chebfun({1, 3, 4, 7, 9}, [-1 -0.5 0 0.25 0.75 1], pref);
g = nextpow2(f);
g_exact = chebfun({0, 2, 2, 3, 4}, [-1 -0.5 0 0.25 0.75 1], pref);
pass(1) = norm(feval(g, x) - feval(g_exact, x), Inf) < ...
    vscale(f)*eps;

f = chebfun(@(x) 4 + 3*sin(x), [-1 0 1], pref);
g = nextpow2(f, pref);
g_exact = @(x) mynextpow2(4 + 3*sin(x));
pass(2) = norm(feval(g, x) - g_exact(x), Inf) < 10*vscale(f)*eps;

% Test for array-valued chebfuns.
f = chebfun([1 3 4 7 9], [-1 1], pref);
g = nextpow2(f);
g_exact = chebfun([0 2 2 3 4], [-1 1], pref);
err = feval(g, x) - feval(g_exact, x);
pass(3) = norm(err(:), Inf) < 10*vscale(f)*eps;

f_op = @(x) 4 + 3*[sin(x) sin(x)];
f = chebfun(@(x) 4 + 3*[sin(x) sin(x)], [-1 0 1], pref);
g = nextpow2(f, pref);
g_exact = @(x) mynextpow2(f_op(x));
err = feval(g, x) - g_exact(x);
pass(4) = norm(err(:), Inf) < 10*vscale(f)*eps;

end

function y = nextpow2Vectorized(x)
%NEXTPOW2VECTORIZED   A manually vectorized NEXTPOW2 for compabibility with
%   versions of MATLAB prior to R2010a.

y = zeros(size(x));
for n = 1:1:numel(x)
    y(n) = nextpow2(x(n));
end

end

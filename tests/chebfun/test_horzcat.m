% Test file for @chebfun/horzcat.m.

function pass = test_horzcat(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

f = chebfun(@sin, [-1 0 1]);
g = chebfun(@cos, [-1 0.5 1]);
h = chebfun(@exp, [-1 -0.5 0 0.5 1]);

Q = [1 f];
Q_exact = @(x) [ones(size(x)) sin(x)];
err = feval(Q, xr) - Q_exact(xr);
pass(1) = numel(f) == 1 && norm(err(:), inf) < 10*vscale(Q)*epslevel(Q);

Q = [f f];
Q_exact = @(x) [sin(x) sin(x)];
err = feval(Q, xr) - Q_exact(xr);
pass(2) = numel(f) == 1 && norm(err(:), inf) < 10*vscale(Q)*epslevel(Q);

Q = [f g];
Q_exact = @(x) [sin(x) cos(x)];
err = feval(Q, xr) - Q_exact(xr);
pass(3) = numel(Q) == 2 && norm(err(:), inf) < 10*vscale(Q)*epslevel(Q);

Q = [f g h];
Q_exact = @(x) [sin(x) cos(x) exp(x)];
err = feval(Q, xr) - Q_exact(xr);
pass(4) = numel(Q) == 3 && norm(err(:), inf) < 10*vscale(Q)*epslevel(Q);

% These functions have the same breakpoints, so we should get an array-valued
% CHEBFUN instead of a quasimatrix.
Q = [f f];
Q = [Q f];
Q_exact = @(x) [sin(x) sin(x) sin(x)];
err = feval(Q, xr) - Q_exact(xr);
pass(5) = numel(Q) == 1 && norm(err(:), inf) < 10*vscale(Q)*epslevel(Q);

end

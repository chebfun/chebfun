% Test file for @chebfun/repmat.m.

function pass = test_repmat(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

f = chebfun(@sin, [-1 0 1]);

Q = repmat(f, 1, 3);
Q_exact = @(x) [sin(x) sin(x) sin(x)];
err = feval(Q, xr) - Q_exact(xr);
pass(1) = norm(err(:), inf) < 10*vscale(Q)*epslevel(Q);

Q = repmat(f, [1 3]);
Q_exact = @(x) [sin(x) sin(x) sin(x)];
err = feval(Q, xr) - Q_exact(xr);
pass(2) = norm(err(:), inf) < 10*vscale(Q)*epslevel(Q);

% TODO:  Comment this in once vertcat() becomes available.
% Q = repmat(f.', 3, 1);
% Q_exact = @(x) [sin(x) sin(x) sin(x)].';
% err = feval(Q, xr) - Q_exact(xr);
% pass(3) = norm(err(:), inf) < 10*vscale(Q)*epslevel(Q);

end

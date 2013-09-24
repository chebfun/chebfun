% Test file for @chebfun/fliplr.m.

function pass = test_flipud(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfun.pref;
end

pref.chebfun.splitting = 1;

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(50, 1) - 1;

% Check the empty input case.
pass(1) = isempty(fliplr(chebfun()));

% Check behavior for row chebfuns.
f = chebfun(@(x) sin(x).*abs(x - 0.1), [-1 0.1 1], pref);
ff = fliplr(f.');
ff_exact = @(x) (sin(-x).*abs(-x - 0.1)).';
err = feval(ff, xr) - ff_exact(xr);
pass(2) = norm(err(:), inf) < 10*ff.vscale.*ff.epslevel;
pass(3) = isequal(fliplr(ff), f.');

g = chebfun(@(x) [sin(x).*abs(x - 0.1) cos(x).*sign(x + 0.2)], ...
    [-1 -0.2 0.1 1], pref);
gg = fliplr(g.');
gg_exact = @(x) [sin(-x).*abs(-x - 0.1) cos(-x).*sign(-x + 0.2)].';
err = feval(gg, xr) - gg_exact(xr);
pass(4) = norm(err(:), inf) < 10*gg.vscale.*gg.epslevel;
pass(5) = isequal(fliplr(gg), g.');

%% Check behavior for column chebfuns.
ff = fliplr(f);
ff_exact = @(x) sin(x).*abs(x - 0.1);
err = feval(ff, xr) - ff_exact(xr);
pass(6) = norm(err(:), inf) < 10*ff.vscale.*ff.epslevel;
pass(7) = isequal(fliplr(ff), f);

gg = fliplr(g);
gg_exact = @(x) [cos(x).*sign(x + 0.2) sin(x).*abs(x - 0.1)];
err = feval(gg, xr) - gg_exact(xr);
pass(8) = norm(err(:), inf) < 10*gg.vscale.*gg.epslevel;
pass(9) = isequal(fliplr(gg), g);

end

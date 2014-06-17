% Test file for @chebfun/flipud.m.

function pass = test_flipud(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

pref.splitting = 1;

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(50, 1) - 1;

% Check the empty input case.
pass(1) = isempty(flipud(chebfun()));

% Check behavior for column chebfuns.
f = chebfun(@(x) sin(x).*abs(x - 0.1), [-1 0.1 1], pref);
ff = flipud(f);
ff_exact = @(x) sin(-x).*abs(-x - 0.1);
err = feval(ff, xr) - ff_exact(xr);
pass(2) = norm(err(:), inf) < 10*ff.vscale.*ff.epslevel;
pass(3) = isequal(flipud(ff), f);

g = chebfun(@(x) [sin(x).*abs(x - 0.1) cos(x).*sign(x + 0.2)], ...
    [-1 -0.2 0.1 1], pref);
gg = flipud(g);
gg_exact = @(x) [sin(-x).*abs(-x - 0.1) cos(-x).*sign(-x + 0.2)];
err = feval(gg, xr) - gg_exact(xr);
pass(4) = norm(err(:), inf) < 10*gg.vscale.*gg.epslevel;
pass(5) = isequal(flipud(gg), g);

% Check behavior for row chebfuns.
ff = flipud(f.');
ff_exact = @(x) (sin(x).*abs(x - 0.1)).';
err = feval(ff, xr) - ff_exact(xr);
pass(6) = norm(err(:), inf) < 10*ff.vscale.*ff.epslevel;
pass(7) = isequal(flipud(ff), f.');

gg = flipud(g.');
gg_exact = @(x) [cos(x).*sign(x + 0.2) sin(x).*abs(x - 0.1)].';
err = feval(gg, xr) - gg_exact(xr);
pass(8) = norm(err(:), inf) < 10*gg.vscale.*gg.epslevel;
pass(9) = isequal(flipud(gg), g.');

%% Test on integration of SINGFUN - A column CHEBFUN:

% Set the domain:
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

pow = -0.5;
op = @(x) (x - dom(1)).^pow.*sin(100*x).*(x - dom(2)).^pow;
f = chebfun(op, dom, 'exps', [pow pow], 'splitting', 'on');
f = f.';
g = flipud(f);
vals_f = feval(f,x);
vals_g = feval(g,x);
err = vals_f - vals_g;
pass(10) = ( norm(err, inf) < get(f,'epslevel')*norm(vals_f, inf) );

%% Test on singular function - a row CHEBFUN:

% Set the domain:
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

pow1 = -0.5;
pow2 = -0.5;
op = @(x) (x - dom(1)).^pow1.*sin(100*x).*(x - dom(2)).^pow2;
opflip = @(x) (dom(1)-x).^pow1.*sin(100*(sum(dom)-x)).*(dom(2)-x).^pow2;

f = chebfun(op, dom, 'exps', [pow1 pow2], 'splitting', 'off');
g = flipud(f);
vals_g = feval(g,x);
vals_exact = feval(opflip,x);
err = vals_g - vals_exact;
pass(11) = ( norm(err, inf) < 1e2*get(f,'epslevel')*norm(vals_exact, inf) );

end

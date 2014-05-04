% Test file for @chebfun/min.m.

function pass = test_min(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
xr = 2 * rand(100, 1) - 1;

% Check empty case.
[y, x] = min(chebfun());
pass(1) = isempty(y) && isempty(x);

% Check operation without breakpoints.
f = chebfun(@(x) ((x - 0.2).^3 - (x - 0.2) + 1).*sec(x - 0.2), pref);
[y, x] = min(f);
y_exact = 0.710869767377087;
pass(2) = all(abs(y - y_exact) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(feval(f, x) - y_exact) <= 10*vscale(f)*epslevel(f));

% Check operation with breakpoints.
f = chebfun(@(x) ((x - 0.2).^3 - (x - 0.2) + 1).*sec(x - 0.2), ...
    linspace(-1, 1, 10), pref);
[y, x] = min(f);
y_exact = 0.710869767377087;
pass(3) = all(abs(y - y_exact) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(feval(f, x) - y_exact) <= 10*vscale(f)*epslevel(f));

%% Check operation for complex-valued chebfuns.
f = chebfun({@(x) exp((1 + 1i)*x), @(x) sec(1i*(x - 0.5))}, [-1 0 1], ...
    pref);
[y, x] = min(f);
y_exact = exp(-1 - 1i);
pass(4) = all(abs(y - y_exact) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(feval(f, x) - y_exact) <= 10*vscale(f)*epslevel(f));

% Check operation for pointValues.
f = chebfun({-1, 1, 2}, [-1, 0, 1, 2]);
[y, ignored] = min(f);
pass(5) = y == -1;

f.pointValues(1,1) = 10;
f.pointValues(3,1) = -10;
[y, x] = min(f);
pass(6) = (y == -10) && (x == 1);

% Check computation of local extrema.
f = chebfun(@(x) sin(x).^2 + sin(x.^2), [0, 4]);
[y, x] = min(f, 'local');
y_exact = [        0
  -0.342247088203205
  -0.971179645473729
   0.284846700239241];
x_exact = [        0
   2.220599667639221
   3.308480466603983
   4.000000000000000];
pass(7) = numel(y == 4) && (norm(y - y_exact, inf) < 10*vscale(f)*epslevel(f));

% Check operation for array-valued chebfuns.
f = chebfun(@(x) [sin(10*x) cos(10*x) exp(x)], [-1 -0.5 0.5 1]);
[y, x] = min(f);
y_exact = [-1 -1 exp(-1)];
fx = feval(f, x(:));
fx = [fx(1, 1) fx(2, 2) fx(3, 3)];
pass(8) = all(abs(y(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(fx(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f));

op = @(x) sin(x).^2 + sin(x.^2);
f = chebfun(@(x) [op(x) op(x/2)], [0, 4]);
[y, x] = min(f, 'local');
y_exact = [        0                  0
  -0.342247088203205  0.070019315123878
  -0.971179645473729  NaN
   0.284846700239241  NaN];
x_exact = [        0                  0
   2.220599667639221  4.000000000000000
   3.308480466603983  NaN
   4.000000000000000  NaN];
fx1 = feval(f, x_exact(:,1));
fx2 = feval(f, x_exact(1:2,2));
pass(9) = isequal(size(y), [4 2]) && ...
    all(isnan(y(3:end,2))) && all(isnan(x(3:end,2)));
pass(10) = norm(y(:,1) - y_exact(:,1), inf) < 10*vscale(f)*epslevel(f) && ...
    norm(y(1:2,2) - y_exact(1:2,2), inf) < 10*vscale(f)*epslevel(f) && ...
    norm(fx1(:,1) - y_exact(:,1), inf) < 10*vscale(f)*epslevel(f) && ...
    norm(fx2(:,2) - y_exact(1:2,2), inf) < 10*vscale(f)*epslevel(f);

% Test min(f, g), where f and g are chebfuns.
f = chebfun(@(x) sin(2*pi*x), [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) cos(2*pi*x), [-1 -0.5 0 0.5 1], pref);
h = min(f, g);
h_exact = @(x) min(sin(2*pi*x), cos(2*pi*x));
pass(11) = norm(feval(h, xr) - h_exact(xr), inf) < 10*vscale(h)*epslevel(h);

g = chebfun(@(x) exp(2*pi*1i*x), [-1 -0.5 0 0.5 1], pref);
h = min(f, g);
h_exact = @(x) min(sin(2*pi*x), exp(2*pi*1i*x));
pass(12) = norm(feval(h, xr) - h_exact(xr), inf) < 10*vscale(h)*epslevel(h);

% NB:  The call to complex() in this next test is to force MATLAB to do
% complex-valued comparison where it wants to do real-valued.  This is
% necessary because g is a complex chebfun, even though one of its columns is
% real.
f = chebfun(@(x) [sin(2*pi*x) cos(2*pi*x)], [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) [exp(2*pi*1i*x) sin(2*pi*x)], [-1 -0.5 0 0.5 1], pref);
h = min(f, g);
h_exact = @(x) [min(sin(2*pi*x), exp(2*pi*1i*x)) ...
    min(cos(2*pi*x), complex(sin(2*pi*x)))];
err = feval(h, xr) - h_exact(xr);
pass(13) = norm(err(:), inf) < 10*vscale(h)*epslevel(h);

% Check 'global' syntax.
f = chebfun(@(x) (x - 0.1).^2 - 1, [-1 -0.5 0 0.5 1], pref);
[y, x] = min(f, 'global');
pass(14) = (abs(y + 1) < 10*epslevel(f)*vscale(f)) && ...
    (abs(feval(f, x) + 1) < 10*epslevel(f)*vscale(f));

% Check error condition.
try
    y = max(f, 'bad');
    pass(15) = false;
catch ME
    pass(15) = strcmp(ME.identifier, 'CHEBFUN:max:flag');
end

%% Check min of a CHEBFUN and a scalar:
f = chebfun(@(x) [sin(x) cos(x)]);
h = min(f, .75);
pass(16) = norm(h([-.9 0 .9].') - [sin(-.9) cos(-.9) ; 0 .75 ; .75 cos(.9)]) ...
    < epslevel(h)*vscale(h);

%% test on function defined on unbounded domain:

% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];

% A blow-up function:
op = @(x) x.*(5+exp(x.^3))./(dom(2)-x);
pref.singPrefs.exponents = [0 -1];
f = chebfun(op, dom, pref); 
[y, x] = min(f);
yExact = -Inf;
xExact = dom(2);
errX = x - xExact;
pass(17) = ( norm(errX, inf) < epslevel(f)*vscale(f) ) && ...
    ( y == yExact );

end

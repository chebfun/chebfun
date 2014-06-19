% Test file for @chebfun/max.m.

function pass = test_max(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
xr = 2 * rand(100, 1) - 1;

% Check empty case.
[y, x] = max(chebfun());
pass(1) = isempty(y) && isempty(x);

% Check operation without breakpoints.
f = chebfun(@(x) ((x - 0.2).^3 - (x - 0.2) + 1).*sec(x - 0.2), pref);
[y, x] = max(f);
y_exact = 1.884217141925336;
pass(2) = all(abs(y - y_exact) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(feval(f, x) - y_exact) <= 10*vscale(f)*epslevel(f));

% Check operation with breakpoints.
f = chebfun(@(x) ((x - 0.2).^3 - (x - 0.2) + 1).*sec(x - 0.2), ...
    linspace(-1, 1, 10), pref);
[y, x] = max(f);
y_exact = 1.884217141925336;
pass(3) = all(abs(y - y_exact) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(feval(f, x) - y_exact) <= 10*vscale(f)*epslevel(f));

% Check operation for complex-valued chebfuns.
f = chebfun({@(x) exp((1 + 1i)*x), @(x) 1-x/10}, [-1 0 1],  pref);
[y, x] = max(f);
y_exact = 1;
pass(4) = all(abs(y - y_exact) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(feval(f, x) - y_exact) <= 10*vscale(f)*epslevel(f));

% Check operation for impulses.
f = chebfun({-1, 1, 2}, [-1, 0, 1, 2]);
[y, ignored] = max(f);
pass(5) = y == 2;

f.pointValues(1,1) = 10;
f.pointValues(3,1) = -10;
[y, x] = max(f);
pass(6) = (y == 10) && (x == -1);

% Check computation of local extrema.
f = chebfun(@(x) sin(x).^2 + sin(x.^2), [0, 4]);
[y, x] = max(f, 'local');
y_exact = [1.923771282655145
           1.117294907913736
           1.343997479566445];
x_exact = [1.323339426259694
           2.781195946808315
           3.776766383330969];
pass(7) = numel(y == 3) && norm(y - y_exact, inf) < 10*vscale(f)*epslevel(f);

% Check operation for array-valued chebfuns.
f = chebfun(@(x) [sin(10*x) cos(10*x) exp(x)], [-1 -0.5 0.5 1]);
[y, x] = max(f);
y_exact = [1 1 exp(1)];
fx = feval(f, x(:));
fx = [fx(1,1) fx(2,2) fx(3,3)];
pass(8) = all(abs(y(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(fx(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f));

f = chebfun(@(x) [exp((1 + 1i)*x) sec(1i*(x - 0.5))], [-1 0 1], ...
    pref);
[y, x] = max(f);
y_exact = [exp(1 + 1i) 1];
fx = feval(f, x(:));
fx = [fx(1,1) fx(2,2)];
pass(9) = all(abs(y(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(fx(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f));

op = @(x) sin(x).^2 + sin(x.^2);
f = chebfun(@(x) [op(x) op(x/2)], [0, 4]);
[y, x] = max(f, 'local');
y_exact = [1.923771282655145  1.923771282655145
           1.117294907913736  NaN
           1.343997479566445  NaN];
x_exact = [1.323339426259694  2.646678852519388
           2.781195946808315  NaN
           3.776766383330969  NaN];
fx1 = feval(f, x_exact(:,1));
fx2 = feval(f, x_exact(1,2));
pass(10) = isequal(size(y), [3 2]) && ...
    all(isnan(y(2:end,2))) && all(isnan(x(2:end,2)));
pass(11) = norm(y(:,1) - y_exact(:,1), inf) < 10*vscale(f)*epslevel(f) && ...
    norm(y(1,2) - y_exact(1,2), inf) < 10*vscale(f)*epslevel(f) && ...
    norm(fx1(:,1) - y_exact(:,1), inf) < 10*vscale(f)*epslevel(f) && ...
    norm(fx2(:,2) - y_exact(1,2), inf) < 10*vscale(f)*epslevel(f);

% Test max(f, g), where f and g are chebfuns.
f = chebfun(@(x) sin(2*pi*x), [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) cos(2*pi*x), [-1 -0.5 0 0.5 1], pref);
h = max(f, g);
h_exact = @(x) max(sin(2*pi*x), cos(2*pi*x));
pass(12) = norm(feval(h, xr) - h_exact(xr), inf) < 10*vscale(h)*epslevel(h);

g = chebfun(@(x) exp(2*pi*1i*x), [-1 -0.5 0 0.5 1], pref);
h = max(f, g);
h_exact = @(x) max(sin(2*pi*x), exp(2*pi*1i*x));
pass(13) = norm(feval(h, xr) - h_exact(xr), inf) < 10*vscale(h)*epslevel(h);

% NB:  The call to complex() in this next test is to force MATLAB to do
% complex-valued comparison where it wants to do real-valued.  This is
% necessary because g is a complex chebfun, even though one of its columns is
% real.
f = chebfun(@(x) [sin(2*pi*x) cos(2*pi*x)], [-1 -0.5 0 0.5 1], pref);
g = chebfun(@(x) [exp(2*pi*1i*x) sin(2*pi*x)], [-1 -0.5 0 0.5 1], pref);
h = max(f, g);
h_exact = @(x) [max(sin(2*pi*x), exp(2*pi*1i*x)) ...
    max(cos(2*pi*x), complex(sin(2*pi*x)))];
err = feval(h, xr) - h_exact(xr);
pass(14) = norm(err(:), inf) < 10*vscale(h)*epslevel(h);

% Check 'global' syntax.
f = chebfun(@(x) -(x - 0.1).^2 + 1, [-1 -0.5 0 0.5 1], pref);
[y, x] = max(f, 'global');
pass(15) = (abs(y - 1) < 10*epslevel(f)*vscale(f)) && ...
    (abs(feval(f, x) - 1) < 10*epslevel(f)*vscale(f));

% Check max(f, [], dim).
f = chebfun(@(x) x);
Q = [f -f];
tol = vscale(Q)*epslevel(Q);

y_a = max(Q, [], 1);
err_a = norm(y_a - [1 1], Inf);
y_q = max(cheb2quasi(Q), [], 1);
err_q = norm(y_q - [1 1], Inf);
y_row = max(Q.', [], 2);
err_row = norm(y_row - [1 1], Inf);
pass(16) = (err_a < tol) && (err_q < tol) && (err_row < tol);

g_a = max(Q, [], 2);
err_a = norm(g_a(xr) - abs(xr), Inf);
g_q = max(cheb2quasi(Q), [], 2);
err_q = norm(g_q(xr) - abs(xr), Inf);
g_row = max(Q.', [], 1);
err_row = norm(g_row(xr) - abs(xr).', Inf);
pass(17) = (err_a < tol) && (err_q < tol) && (err_row < tol);

try
    y = max(Q, [], 3)
    pass(18) = false;
catch ME
    pass(18) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:max:badDim');
end

% Check error condition.
try
    y = max(f, 'bad');
    pass(19) = false;
catch ME
    pass(19) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:max:flag');
end

f = chebfun({@(x) exp((1 + 1i)*x), @(x) sec(1i*(x - 0.5))}, [-1 0 1]);
[y, x] = max(f);
y_exact = 1;
tol = 10*epslevel(f)*vscale(f);
err1 = abs(y - y_exact);
err2 = abs(feval(f, x) - y_exact);
pass(20) = err1 < tol && err2 < tol;

%% Check max of a CHEBFUN and a scalar:
f = chebfun(@(x) [sin(x) cos(x)]);
h = max(f, .75);
pass(21) = norm(h([-.9 0 .8 .9].') - ...
    [.75 .75 ;.75 1 ; .75 .75 ; sin(.9) .75]) < epslevel(h)*vscale(h);

%% Test on function defined on unbounded domain:

% Functions on [a inf]:

% Set the domain:
dom = [1 Inf];

op = @(x) x.*exp(-x);
f = chebfun(op, dom);
[y, x] = max(f);
yExact = exp(-1);
xExact = 1;
errY = y - yExact;
errX = x - xExact;
pass(22) = norm([errY errX], inf) < 100*epslevel(f).*vscale(f);

end

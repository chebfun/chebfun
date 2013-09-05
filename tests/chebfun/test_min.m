% Test file for @chebfun/min.m.

function pass = test_min(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

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

% Check operation for impulses.
f = chebfun({-1, 1, 2}, [-1, 0, 1, 2]);
[y, ignored] = min(f);
pass(5) = y == -1;

f.impulses(1,1,1) = 10;
f.impulses(3,1,1) = -10;
[y, x] = min(f);
pass(6) = (y == -10) && (x == 1);

f = chebfun({-1, 1, 2}, [-1, 0, 1, 2]);
f.impulses(1,1,2) = -1;
f.impulses(2,1,3) = 1;
[y, x] = min(f);
pass(7) = (y == -inf) && (x == -1);

f.impulses(2,1,2) = -1;
[y, x] = min(f);
pass(8) = (y == -inf) && (x == -1);

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
pass(9) = numel(y == 4) && norm(y - y_exact, inf) < 10*vscale(f)*epslevel(f);

% Check operation for array-valued chebfuns.
f = chebfun(@(x) [sin(10*x) cos(10*x) exp(x)], [-1 -0.5 0.5 1]);
[y, x] = min(f);
y_exact = [-1 -1 exp(-1)];
fx = feval(f, x(:));
fx = [fx(1, 1) fx(2, 2) fx(3, 3)];
pass(10) = all(abs(y(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(fx(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f));

f = chebfun({[-1 1 2], [1 -1 3]}, [-1 0 1]);
f.impulses(3, 1, 2) = 1;
f.impulses(2, 3, 3) = -1;
[y, x] = min(f);
pass(11) = isequal(y, [-1 -1 -inf]) && isequal(x, [-1 0 0]);

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
pass(12) = isequal(size(y), [4 2]) && ...
    all(isnan(y(3:end,2))) && all(isnan(x(3:end,2)));
pass(13) = norm(y(:,1) - y_exact(:,1), inf) < 10*vscale(f)*epslevel(f) && ...
    norm(y(1:2,2) - y_exact(1:2,2), inf) < 10*vscale(f)*epslevel(f) && ...
    norm(fx1(:,1) - y_exact(:,1), inf) < 10*vscale(f)*epslevel(f) && ...
    norm(fx2(:,2) - y_exact(1:2,2), inf) < 10*vscale(f)*epslevel(f);

% [TODO]:  Test the min(f, g) syntax, once it is implemented.

end

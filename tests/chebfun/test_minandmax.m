% Test file for @chebfun/minandmax.m.

function pass = test_minandmax(pref)

if ( nargin == 0 )
    pref = chebpref();
end

% Check empty case.
[y, x] = minandmax(chebfun());
pass(1) = isempty(y) && isempty(x);

% Check operation without breakpoints.
f = chebfun(@(x) ((x - 0.2).^3 - (x - 0.2) + 1).*sec(x - 0.2), pref);
[y, x] = minandmax(f);
y_exact = [0.710869767377087 ; 1.884217141925336];
pass(2) = all(abs(y - y_exact) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(feval(f, x) - y_exact) <= 10*vscale(f)*epslevel(f));

% Check operation with breakpoints.
f = chebfun(@(x) ((x - 0.2).^3 - (x - 0.2) + 1).*sec(x - 0.2), ...
    linspace(-1, 1, 10), pref);
[y, x] = minandmax(f);
y_exact = [0.710869767377087 ; 1.884217141925336];
pass(3) = all(abs(y - y_exact) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(feval(f, x) - y_exact) <= 10*vscale(f)*epslevel(f));

% Check operation for complex-valued chebfuns.
f = chebfun({@(x) exp((1 + 1i)*x), @(x) 1-x/10}, [-1 0 1],  pref);
[y, x] = minandmax(f);
y_exact = [exp(-1 - 1i) ; 1];
err1 = abs(y - y_exact);
err2 = abs(feval(f, x) - y_exact);
pass(4) = all(err1 <= 10*vscale(f)*epslevel(f)) && ...
    all(err2 <= 10*vscale(f)*epslevel(f));

% Check operation for impulses.
f = chebfun({-1, 1, 2}, [-1, 0, 1, 2]);
[y, ignored] = minandmax(f);
pass(5) = all(y == [-1 ; 2]);

f.impulses(1,1,1) = 10;
f.impulses(3,1,1) = -10;
[y, x] = minandmax(f);
pass(6) = all(y == [-10 ; 10]) && all(x == [1 ; -1]);

f = chebfun({-1, 1, 2}, [-1, 0, 1, 2]);
f.impulses(1,1,2) = -1;
f.impulses(2,1,3) = 1;
[y, x] = minandmax(f);
pass(7) = all(y == [-inf ; inf]) && all(x == [-1 ; 0]);

f.impulses(2,1,2) = -1;
[y, x] = minandmax(f);
pass(8) = all(y == [-inf ; 2]) && all(x == [-1 ; 1]);

% Check computation of local extrema.
f = chebfun(@(x) sin(x).^2 + sin(x.^2), [0, 4]);
[y, x] = minandmax(f, 'local');
y_exact = [        0
   1.923771282655145
  -0.342247088203205
   1.117294907913736
  -0.971179645473729
   1.343997479566445
   0.284846700239241];
x_exact = [        0
   1.323339426259694
   2.220599667639221
   2.781195946808315
   3.308480466603983
   3.776766383330969
   4.000000000000000];
pass(9) = numel(y == 7) && norm(y - y_exact, inf) < 10*vscale(f)*epslevel(f);

% Check operation for array-valued chebfuns.
f = chebfun(@(x) [sin(10*x) cos(10*x) exp(x)], [-1 -0.5 0.5 1]);
[y, x] = minandmax(f);
y_exact = [-1 -1 exp(-1) ; 1 1 exp(1)];
fx = feval(f, x(:));
fx = [fx(1:2,1) fx(3:4,2) fx(5:6,3)];
pass(10) = all(abs(y(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(fx(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f));

f = chebfun(@(x) [exp((1 + 1i)*x) sec(1i*(x - 0.5))], [-1 0 1], ...
    pref);
[y, x] = minandmax(f);
y_exact = [exp(-1 - 1i) sec(-1.5i) ; exp(1 + 1i) 1];
fx = feval(f, x(:));
fx = [fx(1:2,1) fx(3:4,2)];
pass(11) = all(abs(y(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f)) && ...
    all(abs(fx(:) - y_exact(:)) <= 10*vscale(f)*epslevel(f));

f = chebfun({[-1 1 2], [1 -1 3]}, [-1 0 1]);
f.impulses(3, 1, 2) = 1;
f.impulses(2, 3, 3) = -1;
[y, x] = minandmax(f);
pass(12) = isequal(y, [-1 -1 -inf ; inf 1 3]) && isequal(x, [-1 0 0 ; 1 -1 0]);

op = @(x) sin(x).^2 + sin(x.^2);
f = chebfun(@(x) [op(x) op(x/2)], [0, 4]);
[y, x] = minandmax(f, 'local');
y_exact = [        0                  0
   1.923771282655145  1.923771282655145
  -0.342247088203205  0.070019315123878
   1.117294907913736  NaN
  -0.971179645473729  NaN
   1.343997479566445  NaN
   0.284846700239241  NaN];
x_exact = [        0                  0
   1.323339426259694  2.646678852519388
   2.220599667639221  4.000000000000000
   2.781195946808315  NaN
   3.308480466603983  NaN
   3.776766383330969  NaN
   4.000000000000000  NaN];
fx1 = feval(f, x_exact(:,1));
fx2 = feval(f, x_exact(1:3,2));
pass(13) = isequal(size(y), [7 2]) && ...
    all(isnan(y(4:end,2))) && all(isnan(x(4:end,2)));
pass(14) = norm(y(:,1) - y_exact(:,1), inf) < 10*vscale(f)*epslevel(f) && ...
    norm(y(1:3,2) - y_exact(1:3,2), inf) < 10*vscale(f)*epslevel(f) && ...
    norm(fx1(:,1) - y_exact(:,1), inf) < 10*vscale(f)*epslevel(f) && ...
    norm(fx2(:,2) - y_exact(1:3,2), inf) < 10*vscale(f)*epslevel(f);

%% Test on singular function: piecewise smooth chebfun - splitting on.

dom = [-2 7];
pow = -0.5;
op = @(x) (x - dom(1)).^pow.*(sin(300*x).^2);
pref.singPrefs.exponents = [pow 0];
pref.enableBreakpointDetection = 1;
f = chebfun(op, dom, pref);
[y, x] = minandmax(f);
y_exact = [0 ; Inf];
fx = op(x);
pass(15) = ((max(abs(y - y_exact)) < 2e1*get(f, 'epslevel')) && ... 
          (max(abs(fx - y_exact)) < 2e1*get(f, 'epslevel')));

end

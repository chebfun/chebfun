function pass = test_minandmax(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

f = chebfun({-1,1}, [-1, 0, 1]);
[y, ignored] = minandmax(f);
pass(1) = all(y == [-1 ; 1]);

f = chebfun({-1,1,2}, [-1, 0, 1, 2]);
[y, ignored] = minandmax(f);
pass(2) = all(y == [-1 ; 2]);

f.impulses(1,1,1) = 10;
f.impulses(3,1,1) = -10;
[y, x] = minandmax(f);
pass(3) = all(y == [-10 ; 10]) && all(x == [1 ; -1]);

f = chebfun({-1,1,2}, [-1, 0, 1, 2]);
f.impulses(1,1,2) = -1;
f.impulses(2,1,3) = 1;
[y, x] = minandmax(f);
pass(4) = all(y == [-inf ; inf]) && all(x == [-1 ; 0]);

f.impulses(2,1,2) = -1;
[y, x] = minandmax(f);
pass(5) = all(y == [-inf ; 2]) && all(x == [-1 ; 1]);

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
pass(6) = numel(y == 7) && norm(y - y_exact, inf) < 10*epslevel(f);

% [TODO]: Test array-valued inputs.

    
end
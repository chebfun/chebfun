function pass = test_chebpts3( pref ) 
% Test CHEBPTS3

if ( nargin == 0 ) 
    pref = chebfunpref;
end

tol = 10*pref.cheb3Prefs.chebfun3eps;

% One argument: 
n = 10; 
[xx1, yy1, zz1] = chebfun3.chebpts3(n); 
d = [-1, 1];
x = chebpts(n, d); 
y = chebpts(n, d); 
z = chebpts(n, d);
[xx2, yy2, zz2] = ndgrid(x, y, z); 
pass(1) = norm(xx1(:) - xx2(:)) < tol; 
pass(2) = norm(yy1(:) - yy2(:)) < tol;
pass(3) = norm(zz1(:) - zz2(:)) < tol;

% Three arguments:
m = 10; 
n = 7; 
p = 12;
[xx1, yy1, zz1] = chebfun3.chebpts3(m, n, p);
d = [-1, 1];
x = chebpts(m, d); 
y = chebpts(n, d);
z = chebpts(p, d);
[xx2, yy2, zz2] = ndgrid(x, y, z);
pass(4) = norm(xx1(:) - xx2(:)) < tol; 
pass(5) = norm(yy1(:) - yy2(:)) < tol;
pass(6) = norm(zz1(:) - zz2(:)) < tol;

% Four arguments:
dom = [-2 1 -1 2 -5 -3];
[xx1, yy1, zz1] = chebfun3.chebpts3(m, n, p, dom);
x = chebpts(m, dom(1:2));
y = chebpts(n, dom(3:4));
z = chebpts(p, dom(5:6));
[xx2, yy2, zz2] = ndgrid(x, y, z);
pass(7) = norm(xx1(:) - xx2(:)) < tol; 
pass(8) = norm(yy1(:) - yy2(:)) < tol;
pass(9) = norm(zz1(:) - zz2(:)) < tol;

% Five arguments:
dom = [-2 1 -1 2 -5 -3];
kind = 1;
[xx1, yy1, zz1] = chebfun3.chebpts3(m, n, p, dom, kind);
x = chebpts(m, dom(1:2), kind);
y = chebpts(n, dom(3:4), kind);
z = chebpts(p, dom(5:6), kind);
[xx2, yy2, zz2] = ndgrid(x, y, z);
pass(10) = norm(xx1(:) - xx2(:)) < tol; 
pass(11) = norm(yy1(:) - yy2(:)) < tol;
pass(12) = norm(zz1(:) - zz2(:)) < tol;

end
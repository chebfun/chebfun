function pass = test_prolong(pref)

if ( nargin < 1 )
    pref = funcheb.pref;
end
tol = 10*pref.funcheb.eps;

F = @sin;
f = funcheb1(@(x) F(x), pref);

n = 11;
g = prolong(f, n);
x = funcheb1.chebpts(n);
pass(1) = norm(g.values - F(x), inf) < tol;

g = prolong(f,1);
pass(2) = size(g,1) == 1 && norm(g.values, inf) < tol;

g = prolong(f,2);
y = funcheb1.chebpts(2);
pass(3) = size(g,1) == 2 && norm(g.values - sin(y(1))*[1 ; -1], inf) < tol;

F = @sin;
f = funcheb1(@(x) [F(x), -F(x)], pref);
n = 11;
g = prolong(f, n);
x = funcheb1.chebpts(n);
pass(4) = norm(g.values - [F(x), -F(x)], inf) < tol;

g = prolong(f,1);
pass(5) = size(g,1) == 1 && norm(g.values, inf) < tol;

g = prolong(f,2);
pass(6) = size(g,1) == 2 && norm(g.values - sin(y(1))*[1 -1 ; -1 1], inf) < tol;

g = prolong(f, length(f));
pass(7) = all(f.values(:) == g.values(:));

n = 32;
g = prolong(f, n);
x = funcheb1.chebpts(n);
pass(8) = norm(g.values - [F(x), -F(x)], inf) < tol;

n = 100;
g = prolong(f, n);
x = funcheb1.chebpts(n);
pass(9) = norm(g.values - [F(x), -F(x)], inf) < tol;

F = @(x) sin(1000*x);
f = funcheb1(@(x) [F(x), -F(x)], pref);
n = 32;
g = prolong(f, n);
x = funcheb1.chebpts(n);
pass(10) = norm(g.values - [F(x), -F(x)], inf) < length(f)*tol;

n = 100;
g = prolong(f, n);
x = funcheb1.chebpts(n);
pass(11) = norm(g.values - [F(x), -F(x)], inf) < length(f)*tol;

g = prolong(f,1);
pass(12) = size(g, 1) == 1 && norm(g.values, inf) < length(f)*tol;

g = prolong(f,2);
pass(13) = size(g, 1) == 2 && norm(g.values - sin(1000*y(1))*[1 -1 ; -1 1], inf) < length(f)*tol;

F = @(x) cos(1000*x);
f = funcheb1(@(x) [F(x), -F(x)], pref);

g = prolong(f,1);
pass(14) = size(g, 1) == 1 && norm(g.values - [1, -1], inf) < tol;

g = prolong(f,2);
pass(15) = length(g) == 2 && norm(g.values - cos(1000*y(1))*[1 -1 ; 1 -1], inf) < length(f)*tol;

v = [1 2 3];
f = funcheb1(v, pref);
g = prolong(f, 5);
pass(16) = norm(g.values - repmat([1 2 3], 5, 1), inf) < tol;

end

function pass = test_prolong(pref)

if ( nargin < 1 )
    pref = funcheb2.pref;
end
tol = 10*pref.funcheb2.eps;

F = @sin;
f = funcheb2(@(x) F(x), [], pref);

n = 11;
g = prolong(f, n);
x = funcheb2.chebpts(n);
pass(1) = norm(g.values - F(x), inf) < tol;

g = prolong(f,1);
pass(2) = size(g,1) == 1 && norm(g.values, inf) < tol;

g = prolong(f,2);
pass(3) = size(g,1) == 2 && norm(g.values - sin(1)*[-1 ; 1], inf) < tol;

F = @sin;
f = funcheb2(@(x) [F(x), -F(x)], [], pref);
n = 11;
g = prolong(f, n);
x = funcheb2.chebpts(n);
pass(4) = norm(g.values - [F(x), -F(x)], inf) < tol;

g = prolong(f,1);
pass(5) = size(g,1) == 1 && norm(g.values, inf) < tol;

g = prolong(f,2);
pass(6) = size(g,1) == 2 && norm(g.values - sin(1)*[-1 1 ; 1 -1], inf) < tol;

g = prolong(f, length(f));
pass(7) = all(f.values(:) == g.values(:));

n = 32;
g = prolong(f, n);
x = funcheb2.chebpts(n);
pass(8) = norm(g.values - [F(x), -F(x)], inf) < tol;

n = 100;
g = prolong(f, n);
x = funcheb2.chebpts(n);
pass(9) = norm(g.values - [F(x), -F(x)], inf) < tol;

F = @(x) sin(1000*x);
f = funcheb2(@(x) [F(x), -F(x)], [], pref);
n = 32;
g = prolong(f, n);
x = funcheb2.chebpts(n);
pass(10) = norm(g.values - [F(x), -F(x)], inf) < length(f)*tol;

n = 100;
g = prolong(f, n);
x = funcheb2.chebpts(n);
pass(11) = norm(g.values - [F(x), -F(x)], inf) < length(f)*tol;

g = prolong(f,1);
pass(12) = size(g, 1) == 1 && norm(g.values, inf) < tol;

g = prolong(f,2);
pass(13) = size(g, 1) == 2 && norm(g.values - sin(1000)*[-1 1 ; 1 -1], inf) < tol;

F = @(x) cos(1000*x);
f = funcheb2(@(x) [F(x), -F(x)], [], pref);

g = prolong(f,1);
pass(14) = size(g, 1) == 1 && norm(g.values - [1, -1], inf) < tol;

g = prolong(f,2);
pass(15) = length(g) == 2 && norm(g.values - cos(1000)*[1 -1 ; 1 -1], inf) < tol;

end

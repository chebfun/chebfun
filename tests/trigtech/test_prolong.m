% Test file for trigtech/prolong.m

function pass = test_prolong(pref)

if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

F = @(x) exp(sin(pi*x));
f = testclass.make(F, [], pref);

% Odd prolongation.
k = 101;
go = prolong(f, k);
x = testclass.trigpts(k);
pass(1) = size(go,1) == k  && norm(go.values - F(x), inf) < 10*go.vscale.*go.epslevel;

% Even prolongation.
k = 100;
ge = prolong(f, k);
x = testclass.trigpts(k);
pass(2) = size(ge,1) == k  && norm(ge.values - F(x), inf) < 10*ge.vscale.*ge.epslevel;

% Odd restriction from odd length f.
k = 89;
g = prolong(go, k);
x = testclass.trigpts(k);
pass(3) = size(g,1) == k  && norm(g.values - F(x), inf) < 10*g.vscale.*g.epslevel;

% Even restriction from odd length f.
k = 88;
g = prolong(go, k);
x = testclass.trigpts(k);
pass(4) = size(g,1) == k  && norm(g.values - F(x), inf) < 10*g.vscale.*g.epslevel;

% Odd restriction from even length f.
k = 89;
g = prolong(ge, k);
x = testclass.trigpts(k);
pass(5) = size(g,1) == k  && norm(g.values - F(x), inf) < 10*g.vscale.*g.epslevel;

% Even restriction from even length f.
k = 88;
g = prolong(ge, k);
x = testclass.trigpts(k);
pass(6) = size(g,1) == k  && norm(g.values - F(x), inf) < 10*g.vscale.*g.epslevel;

% Restriction to length 1.
g = prolong(f, 1);
pass(7) = size(g,1) == 1;

% Array-valued tests.
F = @(x) exp(sin(pi*x));
f = testclass.make(@(x) [F(x), -F(x)], [], pref);
k = 101;
g = prolong(f, k);
x = testclass.trigpts(k);
values = g.coeffs2vals(g.coeffs);
pass(8) = size(g,1) == k && norm(values - [F(x), -F(x)], inf) < ...
    10*max(g.vscale.*g.epslevel);

g = prolong(f, 1);
pass(9) = size(g,1) == 1;

% Same length test
g = prolong(f, size(f,1));
fvalues = f.coeffs2vals(f.coeffs);
gvalues = g.coeffs2vals(g.coeffs);
pass(10) = all(fvalues(:) == gvalues(:));

% Values test.
v = [1 2 3];
f = testclass.make(v, [], pref);
g = prolong(f, 5);
values = g.coeffs2vals(g.coeffs);
pass(11) = norm(values - repmat([1 2 3], 5, 1), inf) < ...
   10*max(g.vscale.*g.epslevel);

end

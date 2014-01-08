% Test file for @chebfun/sum.m.

function pass = test_sum(pref)

% Obtain preferences.
if ( nargin == 0 )
    pref = chebpref();
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% Check the empty case.
pass(1) = sum(chebfun()) == 0;

% Check operation in the general case.
f = chebfun({@(x) exp(4*pi*1i*x), @exp, @exp}, [-1 0 0.5 1], pref);
pass(2) = abs(sum(f) - (exp(1) - 1)) < 10*vscale(f)*epslevel(f);

% Check operation for row chebfuns.
ft = f.';
pass(3) = abs(sum(ft) - (exp(1) - 1)) < 10*vscale(ft)*epslevel(ft);

% Check operation with impulses.
f.impulses(3,1,2) = 2;
pass(4) = abs(sum(f) - (exp(1) + 1)) < 10*vscale(f)*epslevel(f);

% Check sum over a subdomain.
pass(5) = abs(sum(f, [-1 1]) - (exp(1) + 1)) < 10*vscale(f)*epslevel(f);
pass(6) = abs(sum(f, [-1 0])) < 10*vscale(f)*epslevel(f);
pass(7) = abs(sum(f, [0 1]) - (exp(1) + 1)) < 10*vscale(f)*epslevel(f);

% Check sum between chebfun limits.
f = chebfun(@exp, [-1 -0.5 0 0.5 1], pref);
a = chebfun(@(x) x.^2 - 1, [-1 1]);
b = chebfun(@(x) -x.^2 + 1, [-1 1]);

F1 = sum(f, a, 1);
F1_exact = @(x) exp(1) - exp(x.^2 - 1);
pass(8) = norm(feval(F1, xr) - F1_exact(xr), inf) < 10*vscale(F1)*epslevel(F1);

F2 = sum(f, -1, b);
F2_exact = @(x) exp(-x.^2 + 1) - exp(-1);
pass(9) = norm(feval(F2, xr) - F2_exact(xr), inf) < 10*vscale(F2)*epslevel(F2);

F3 = sum(f, a, b);
F3_exact = @(x) exp(-x.^2 + 1) - exp(x.^2 - 1);
pass(10) = norm(feval(F3, xr) - F3_exact(xr), inf) < 10*vscale(F3)*epslevel(F3);

% Check operation for array-valued chebfuns.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1]);
pass(11) = norm(sum(f) - [0 2*sin(1) (exp(1) - exp(-1))], inf) < ...
    10*vscale(f)*epslevel(f);

ft = f.';
pass(12) = norm(sum(ft) - [0 2*sin(1) (exp(1) - exp(-1))], inf) < ...
    10*vscale(ft)*epslevel(ft);

g = f;
g.impulses(2,2,2) = 1;
pass(13) = norm(sum(g) - [0 (2*sin(1) + 1) (exp(1) - exp(-1))], inf) < ...
    10*vscale(g)*epslevel(g);

pass(14) = norm(sum(f, [-1 1]) - [0 2*sin(1) (exp(1) - exp(-1))], inf) < ...
    10*vscale(f)*epslevel(f);

pass(15) = norm(sum(f, [-1 0]) - [(cos(-1) - 1) sin(1) (1 - exp(-1))], inf) ...
    < 10*vscale(f)*epslevel(f);

F1 = sum(f, a, 1);
F1_col1_exact = @(x) cos(x.^2 - 1) - cos(1);
F1_col2_exact = @(x) sin(1) - sin(x.^2 - 1);
F1_col3_exact = @(x) exp(1) - exp(x.^2 - 1);
F1_exact = @(x) [F1_col1_exact(x) F1_col2_exact(x) F1_col3_exact(x)];
err = feval(F1, xr) - F1_exact(xr);
pass(16) = norm(err(:), inf) < 10*vscale(F1)*epslevel(F1);

F2 = sum(f, -1, b);
F2_col1_exact = @(x) cos(-1) - cos(-x.^2 + 1);
F2_col2_exact = @(x) sin(-x.^2 + 1) - sin(-1);
F2_col3_exact = @(x) exp(-x.^2 + 1) - exp(-1);
F2_exact = @(x) [F2_col1_exact(x) F2_col2_exact(x) F2_col3_exact(x)];
err = feval(F2, xr) - F2_exact(xr);
pass(17) = norm(err(:), inf) < 10*vscale(F2)*epslevel(F2);

F3 = sum(f, a, b);
F3_col1_exact = @(x) -cos(-x.^2 + 1) + cos(x.^2 - 1);
F3_col2_exact = @(x) sin(-x.^2 + 1) - sin(x.^2 - 1);
F3_col3_exact = @(x) exp(-x.^2 + 1) - exp(x.^2 - 1);
F3_exact = @(x) [F3_col1_exact(x) F3_col2_exact(x) F3_col3_exact(x)];
err = feval(F3, xr) - F3_exact(xr);
pass(18) = norm(err(:), inf) < 10*vscale(F3)*epslevel(F3);

% Check dim argument.
g = sum(f, 2);
g_exact = @(x) sin(x) + cos(x) + exp(x);
pass(19) = norm(feval(g, xr) - g_exact(xr), inf) < 10*vscale(g)*epslevel(g);

g = sum(ft, 1);
g_exact = @(x) (sin(x) + cos(x) + exp(x)).';
pass(20) = norm(feval(g, xr) - g_exact(xr), inf) < 10*vscale(g)*epslevel(g);

% Check error conditions.
try
    s = sum(f, -2, 2);
    pass(21) = false;
catch ME
    pass(21) = strcmp(ME.identifier, 'CHEBFUN:sum:ab');
end

try
    s = sum(f, -2, b);
    pass(22) = false;
catch ME
    pass(22) = strcmp(ME.identifier, 'CHEBFUN:sum:a');
end

try
    s = sum(f, a, 2);
    pass(23) = false;
catch ME
    pass(23) = strcmp(ME.identifier, 'CHEBFUN:sum:b');
end

end

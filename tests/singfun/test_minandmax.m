% Test file for singfun/minandmax.m

function pass = test_minandmax(pref)

% Get preferences.
if ( nargin < 1 )
    pref = singfun.pref;
end

% Set a tolerance.
tol = 10*pref.singfun.eps;

% Generate a few random points to use as test values.
seedRNG(6178);
d = 2;
x = 2*(1-10^(-d)) * rand(100, 1) - (1-10^(-d));

% The order of the exponents:
a = 0.56;
b = -0.56;
c = 1.28;
d = -1.28;

% Pre-allocate pass matrix
pass = zeros(1, 6);

% fractional root at the left endpoint and the function value is bounded.
f = singfun(@(x) (1+x).^a.*exp(x), [a 0], {'sing', 'none'}, pref);
[y, x] = minandmax(f);
y_exact = [0; 2^a*exp(1)];
x_exact = [-1; 1];
err_x = norm(x-x_exact, inf);
err_y = norm(y-y_exact, inf);
pass(1) = (max([err_x err_y]) < tol*normest(f));

% fractional pole at the left endpoint and the function value is unbounded.
f = singfun(@(x) (1+x).^d.*sin(50*pi*x), [d+1 0], {'sing', 'none'}, pref);
[y, x] = minandmax(f);
y_exact = [-1 Inf];
err_y = y(1) - y_exact(1);
pass(2) = (abs(err_y) < tol*normest(f)) && (y(2) == Inf);

% fractional root at the right endpoint and the smooth part has no roots in [-1 1].
f = singfun(@(x) (1-x).^c.*cos(x), [0 c], {'none', 'sing'}, pref);
r = roots(f);
r_exact = 1;
err = r - r_exact;
pass(3) = (norm(err, inf) < tol*normest(f));

% no fractional pole but a root at the left endpoint.
f = singfun(@(x) (1-x).^b.*(exp(x)-exp(1)), [0 1+b], {'none', 'sing'}, pref);
r = roots(f);
r_exact = 1;
err = r - r_exact;
pass(4) = (norm(err, inf) < tol*normest(f));

% a combination of fractional pole and fractional root.
f = singfun(@(x) (1+x).^b.*sin(x).*(1-x).^c, [b c], {'sing', 'sing'}, pref);
r = roots(f);
r_exact = [0 1];
err = r - r_exact;
pass(5) = (norm(err, inf) < tol*normest(f));

% Check the case with roots close to endpoints.
p = 1-1e-14;
f = singfun(@(x) (1+x).^b.*sin(x-p).*(1-x).^b, [b b], {'sing', 'sing'}, pref);
r = roots(f);
r_exact = p;
err = r - r_exact;
pass(6) = (norm(err, inf) < tol*normest(f));

end
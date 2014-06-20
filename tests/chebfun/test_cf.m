% Test function for @chebfun/cf.m.
% Based on the Chebfun v4 test written by LNT & Joris Van Deun, Dec. 4, 2009.

function pass = test_cf(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

f = chebfun(@exp, pref);
[p, q, r, lam] = cf(f, 2);
pass(1) = abs(lam - 0.045017) < 1e-4;

f = chebfun(@(x) exp(-1 + x/2), [0 4], pref);
[p, q, r, lam] = cf(f, 2);
pass(2) = abs(lam - 0.045017) < 1e-4;

x = chebfun(@(x) x, [-1 1], pref);

f = cos(x);
[p, q] = cf(f, 1, 1);
pass(3) = abs(p(.3) - .77015046914) < 1e-4;

f = abs(x);
[p, q] = cf(f, 4, 4, 100);
pass(4) = norm(f - p./q) < .05;

f = exp(exp(x));
[p, q, r] = cf(f, 0, 10);
xx = linspace(-1, 1, 17);
pass(5) = norm(f(xx) - r(xx)) < 1e-4;

% Test an example that is not based on the domain [-1, 1].
f = chebfun(@exp, [2 6], pref);
[p, q, r] = cf(f, 5, 5);
xx = linspace(2, 6, 100);
pass(6) = norm(feval(f, xx) - r(xx), inf) < 1e-6;

% Test CF with quasimatrix and array-valued input to ensure the output is of
% the correct form.
f = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 1], pref);
[p, q, r, s] = cf(f, 2, 2);
pass(7) = (numel(p) == 3) && (numel(q) == 3) && (numel(r) == 3) && ...
    (numel(s) == 3);

f = cheb2quasi(f).';
[p, q, r, s] = cf(f, 2, 2);
pass(8) = (numel(p) == 3) && (numel(q) == 3) && (numel(r) == 3) && ...
    (numel(s) == 3) && p(1,:).isTransposed && q(1,:).isTransposed && ...
    iscolumn(r) && iscolumn(s);

end

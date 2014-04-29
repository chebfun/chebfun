% Test function for @chebfun/cf.m.
% Based on the Chebfun v4 test written by LNT & Joris Van Deun, Dec. 4, 2009.

function pass = test_cf(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

f = chebfun(@exp);
[p, q, r, lam] = cf(f, 2);
pass(1) = abs(lam - 0.045017) < 1e-4;

f = chebfun(@(x) exp(-1 + x/2), [0 4]);
[p, q, r, lam] = cf(f, 2);
pass(2) = abs(lam - 0.045017) < 1e-4;

x = chebfun(@(x) x, [-1 1]);

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
f = chebfun(@exp, [2 6]);
[p, q, r] = cf(f, 5, 5);
xx = linspace(2, 6, 100);
pass(6) = norm(feval(f, xx) - r(xx), inf) < 1e-6;

end

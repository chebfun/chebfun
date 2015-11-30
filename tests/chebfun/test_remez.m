% Test file for @chebfun/remez.m.
% Based on the Chebfun v4 test written by LNT, March. 27, 2009.

function pass = test_remez(pref)

% Generate a few random points to use as test values.
seedRNG(6178);
xx = 2 * rand(100, 1) - 1;

x = chebfun(@(x) x, [-1 1]);

% Test known exact answer.
f = abs(x) + x;
pexact = .5 + x;
pbest = remez(f, 1);
err = norm(pbest - pexact);
pass(1) = err < 1e-10;

% Test remez() and cf().
f = exp(sin(exp(x)));
pcf = cf(f, 7);
pbest = remez(f, 7);
pass(2) = norm(pcf - pbest) < 0.0003;

% Test rational mode.
f = ((x + 3).*(x - 0.5))./(x.^2 - 4);
[pbest, qbest, rbest] = remez(f, 2, 2, 'tol', 1e-12, 'maxIter', 20);
pass(3) =  norm(f(xx) - rbest(xx), inf) < 1e-10;

% Test an example with a breakpoint.
x = chebfun(@(x) x, [-1 0 1]);
f = ((x - 3).*(x + 0.2).*(x - 0.7))./((x - 1.5).*(x + 2.1));
[pbest, qbest, rbest] = remez(f, 3, 2);
pass(4) =  norm(f(xx) - rbest(xx), inf) < 1e-10;

% Test an example from ATAP.
x = chebfun(@(x) x);
f = abs(x);
pbest = remez(f, 3);
pbest_exact = x.^2 + 1/8;
pass(5) = norm(pbest(xx) - pbest_exact(xx), inf) < 10*eps;

% Test a zero function (#1656)
f = 0*x;
pbest = remez(f, 2);
pbest_exact = 0*x;
pass(6) = norm(pbest(xx) - pbest_exact(xx), inf) < 10*eps;

end

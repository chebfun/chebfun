% Test file for @chebfun/remez.m.
% Based on the Chebfun v4 test written by LNT, March. 27, 2009.

function pass = remeztest(pref)

x = chebfun(@(x) x, [-1 1]);

% Test known exact answer.
f = abs(x) + x;
pexact = .5 + x;
pbest = remez(f, 1);
err = norm(pbest - pexact);
pass(1) = err < 1e-10;

% Test remez and cf.
f = exp(sin(exp(x)));
pcf = cf(f, 7);
pbest = remez(f, 7);
pass(2) = norm(pcf - pbest) < 0.0003;

end

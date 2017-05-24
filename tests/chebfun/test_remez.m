% Test file for @chebfun/remez.m.

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
[~, ~, rbest] = remez(f, 2, 2, 'tol', 1e-12, 'maxIter', 20);
pass(3) =  norm(f(xx) - rbest(xx), inf) < 1e-10;

% Test an example with a breakpoint.
x = chebfun(@(x) x, [-1 0 1]);
f = ((x - 3).*(x + 0.2).*(x - 0.7))./((x - 1.5).*(x + 2.1));
[~, ~, rbest] = remez(f, 3, 2);
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

% Tests giving trouble on 26 Aug 2016:
f = abs(x);
[p,q] = remez(f,0,0);
err = norm(f-p./q,inf);
pass(7) = (abs(err-.5) < 1e-10);

[p,q] = remez(f,2,2);
err = norm(f-p./q,inf);
pass(8) = (abs(err-.043689) < 1e-3);


% Make sure that it works for m=0 and odd f.
f = x.^3;
try
    [~,~] = remez(f,0,2);
    pass(9) = 1;
catch ME 
    pass(9) = false;
end


% A function with huge amplitude
f1 = chebfun('exp(x)'); p1 = remez(f1,7); err1 = norm(f1-p1,inf);
f2 = chebfun('1e100*exp(x)'); p2 = remez(f2,7); err2 = norm(f2-p2,inf);
pass(10) = abs(err1-err2/1e100) < 1e-3;


% A function with tiny amplitude
f1 = chebfun('exp(x)'); [~,~,r1] = remez(f1,1,3); sample1 = r1(.3)-f1(.3);
f2 = chebfun('1e-100*exp(x)'); [~,~,r2] = remez(f2,1,3); sample2 = r2(.3)-f2(.3);
pass(11) = abs(sample1-sample2*1e100) < 1e-3;


% A function on a big domain
x = chebfun('x',1e30*[-1 2]);
f = abs(x);
p = remez(f,30);
pass(12) = (norm(f-p,inf)/1e30 - .0135210) < .01;

% Check correct output formatting for polynomial Remez
f = chebfun('exp(x)');
[p,err] = remez(f,7);
pass(13) = (abs(norm(f-p,inf)-err) < 1e-14);

end

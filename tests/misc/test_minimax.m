% Test file for minimax.m.

function pass = test_minimax( pref )

% Generate a few random points to use as test values.
seedRNG(6178);
xx = 2 * rand(100, 1) - 1;

x = chebfun(@(x) x, [-1 1]);

% Test known exact answer.
f = abs(x) + x;
pexact = .5 + x;
pbest = minimax(f, 1);
err = norm(pbest - pexact);
pass(1) = err < 1e-10;

% Test minimax() and cf().
f = exp(sin(exp(x)));
pcf = cf(f, 7);
pbest = minimax(f, 7);
pass(2) = norm(pcf - pbest) < 0.0003;

% Test rational mode.
f = ((x + 3).*(x - 0.5))./(x.^2 - 4);
[~, ~, rbest] = minimax(f, 2, 2, 'tol', 1e-12, 'maxIter', 20);
pass(3) =  norm(f(xx) - rbest(xx), inf) < 1e-10;

% Test an example with a breakpoint.
x = chebfun(@(x) x, [-1 0 1]);
f = ((x - 3).*(x + 0.2).*(x - 0.7))./((x - 1.5).*(x + 2.1));
[~, ~, rbest] = minimax(f, 3, 2);
pass(4) =  norm(f(xx) - rbest(xx), inf) < 1e-10;

% Test an example from ATAP.
x = chebfun(@(x) x);
f = abs(x);
pbest = minimax(f, 3);
pbest_exact = x.^2 + 1/8;
pass(5) = norm(pbest(xx) - pbest_exact(xx), inf) < 10*eps;

% Test a zero function (#1656)
f = 0*x;
pbest = minimax(f, 2);
pbest_exact = 0*x;
pass(6) = norm(pbest(xx) - pbest_exact(xx), inf) < 10*eps;

% Tests giving trouble on 26 Aug 2016:
f = abs(x);
[p,q] = minimax(f,0,0);
err = norm(f-p./q,inf);
pass(7) = (abs(err-.5) < 1e-10);

[p,q] = minimax(f,2,2);
err = norm(f-p./q,inf);
pass(8) = (abs(err-.043689) < 1e-3);


% Make sure that it works for m=0 and odd f.
f = x.^3;
try
    [~,~] = minimax(f,0,2);
    pass(9) = 1;
catch ME
    pass(9) = false;
end


% A function with huge amplitude
f1 = chebfun('exp(x)'); p1 = minimax(f1,7); err1 = norm(f1-p1,inf);
f2 = chebfun('1e100*exp(x)'); p2 = minimax(f2,7); err2 = norm(f2-p2,inf);
pass(10) = abs(err1-err2/1e100) < 1e-3;


% A function with tiny amplitude
f1 = chebfun('exp(x)');
[~,~,r1] = minimax(f1,1,3); sample1 = r1(.3)-f1(.3);
f2 = chebfun('1e-100*exp(x)');
[~,~,r2] = minimax(f2,1,3); sample2 = r2(.3)-f2(.3);
pass(11) = abs(sample1-sample2*1e100) < 1e-3;


% A function on a big domain
x = chebfun('x',1e30*[-1 2]);
f = abs(x);
p = minimax(f,30);
pass(12) = (norm(f-p,inf)/1e30 - .0135210) < .01;

% Higher degree abs(x) approximation
x = chebfun('x');
f = abs(x);
[~, ~, ~,err,~] = minimax(f, 30, 30, 'silent');
pass(13) = abs(err-2.1739878e-7)/2.1739878e-7 < 1e-3;

% Check correct output formatting for polynomial Remez
f = chebfun('exp(x)');
[p,err] = minimax(f,7);
pass(14) = (abs(norm(f-p,inf)-err) < 1e-14);

% Test the usage of function handle and string inputs
[p1,err1] = minimax('exp(x)',5);
[p2,err2] = minimax(@(x) exp(x), 5);
[p3,err3] = minimax(chebfun('exp(x)'), 5);
pass(15) = ((err1 == err2) && (abs(err1 - err3) < 1e-15));

% Test the use of minimax on a singfun object (#1405)
x = chebfun('x');
f = sqrt(abs(x-.1));
[p,err] = minimax(f, 5);
xx = linspace(-1,1,10000);
norme1 = max(abs(f(xx)-p(xx)));
pass(16) = (abs(err - norme1)/err < 1e-4);

[~,~,r,err,~] = minimax(f,4,4,'silent');
norme2 = max(abs(f(xx)-r(xx)));
pass(17) = (abs(err - norme2)/err < 1e-4);

x = chebfun('x'); f = 1e40*abs(x);
[p,q,rh,err] = minimax(f,5,5);
pass(18) = (err < 1e38);

% Test poles and zeros of the best approximation

[p,q,~,~,status] = minimax(@(x) sqrt(x), [0,1], 4,4);
zer1 = roots(p,'all'); zer1 = sort(zer1); zer2 = sort(status.zer);
pol1 = roots(q,'all'); pol1 = sort(pol1); pol2 = sort(status.pol);
pass(19) = ( norm(zer1-zer2,Inf)/norm(zer1,Inf) < 1e-5 && ...
    norm(pol1-pol2,Inf)/norm(pol1,Inf) < 1e-5 );

% Test for issue #2244
[p1, err1] = minimax(@exp, 5);
[p2, err2] = minimax(@(x) exp(x), 5);
pass(20) = (err1 == err2);

end

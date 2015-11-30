function pass = test_abs(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebfunpref();
end

% Initialise random vectors:
seedRNG(6178);
x = 2 * rand(100, 1) - 1;
x2 = 101 * rand(300, 1) - 1;
x3 = 4 * rand(100, 1) - 2;

%% Simple tests

% This should not introduce a break at zero:
f = chebfun('x.^2', pref);
f1 = abs(f);
pass(1,1) = (numel(f1.funs) == 1) || (numel(f1.funs) == 2); % Either is OK.

% Test positivity:
f = chebfun(@(x) cos(3*x), pref);
f1 = abs(f);
pass(1,2) = all(feval(f1, x) >= 0);

% Test if pointValues are dealt with correctly: 
f = restrict(f + 2, sort(x)');
%f = chebfun(@(x) sin(x) + 2, sort(x)', pref);
f.pointValues = x;
f1 = abs(f);
pass(1,3) = all(f1.pointValues == abs(x));

% Test also on longer intervals:
f = chebfun(@(x) x, [-1 100], pref);
f1 = abs(f);
pass(1,4) = all(feval(f1, x2) >= 0);

% Test functions with breakpoints:
%f = chebfun(@(x) x, [-1 1e-10 1], pref);  % does not work, but should !?
f = chebfun(@(x) x, [-1 eps 1], pref);
f1 = abs(f);
f = chebfun(@(x) abs(x), [-1 0 1], pref);
tol = 10*eps*get(f, 'hscale');
pass(1,5) = normest(f - f1) < tol;

% The absolute value of all points on the unit circle should be 1:
f = chebfun(@(x) exp(1i*x), pref);
g = abs(f);
tol = 10*eps*get(g, 'hscale');
pass(1,6) = normest(g-1) < tol ;

% Real, imaginary and complex CHEBFUN objects:
fHandle = {@(x) sin(pi*x), @(x) -sin(pi*x), @(x) sin(pi*x), @(x) -sin(pi*x)};
f = chebfun(fHandle, -2:2, pref);
tol = 20*eps*get(f, 'hscale');

%% Real CHEBFUN
gHandle1 = @(x) sin(pi*x);
g = chebfun(@(x) gHandle1(x), -2:2, pref);
g1 = abs(g);
pass(2,1) = length(g1.funs) == 4;
pass(2,2) = normest(f - g1) < tol;
pass(2,3) = all(g1.pointValues == abs(sin(pi*(-2:2)))');
h1 = chebfun(@(x) abs( gHandle1(x) ), -2:2, pref);
pass(2,4) = length(h1.funs) == 4;
pass(2,5) = normest(f - h1) < tol;
pass(2,6) = all(h1.pointValues == abs(sin(pi*(-2:2)))');

%% Imaginary CHEBFUN
gHandle2 = @(x) 1i*cos(pi*(x-.5));
g = chebfun(@(x) gHandle2(x), -2:2, pref);
g2 = abs(g);
pass(3,1) = length(g2.funs) == 4;
pass(3,2) = normest(f - g2) < tol;
pass(3,3) = isreal(g2);
h2 = chebfun(@(x) abs( gHandle2(x) ), -2:2, pref);
pass(3,4) = length(h2.funs) == 4;
pass(3,5) = normest(f - h2) < tol;
pass(3,6) = isreal(h2);

%% Complex CHEBFUN
f1 = sqrt(2)*f;
g = chebfun(@(x) gHandle1(x) + gHandle2(x), -2:2, pref);
g3 = abs(g);
pass(4,1) = length(g3.funs) == 4;
pass(4,2) = normest(f1 - g3) < tol;
pass(4,3) = isreal(g3);
h3 = chebfun(@(x) abs(gHandle1(x) + gHandle2(x)), -2:2, pref);
pass(4,4) = length(h3.funs) == 4;
pass(4,5) = normest(f1 - h3) < tol;
pass(4,6) = isreal(h3);

%% Array-valued CHEBFUN
f1 = chebfun(@(x) feval(f, [x, x]) , -2:2, pref);
gHandle2 = @(x) cos(pi*(x-.5));
g = chebfun(@(x) [gHandle1(x), gHandle2(x)] , -2:2, pref);
g4 = abs(g);
pass(5,1) = length(g4.funs) == 4;
pass(5,2) = normest(f1 - g4) < tol;
pass(5,3) = all(all(feval(g4, x3) >= 0));
h4 = chebfun(@(x) abs([gHandle1(x), gHandle2(x)]), -2:2, pref);
pass(5,4) = length(h4.funs) == 4;
pass(5,5) = normest(f1 - h4) < tol;
pass(5,6) = all(all(feval(h4, x3) >= 0));

%% CHEBFUN Quasimatrix
f1 = chebfun(@(x) feval(f, [x, x]) , -2:2, pref);
f1 = quasi2cheb(f1);
gHandle2 = @(x) cos(pi*(x-.5));
g = chebfun(@(x) [gHandle1(x), gHandle2(x)] , -2:2, pref);
g = quasi2cheb(g);
g4 = abs(g);
pass(6,1) = length(g4.funs) == 4;
pass(6,2) = normest(f1 - g4) < tol;
pass(6,3) = all(all(feval(g4, x3) >= 0));
h4 = chebfun(@(x) abs([gHandle1(x), gHandle2(x)]), -2:2, pref);
pass(6,4) = length(h4.funs) == 4;
pass(6,5) = normest(f1 - h4) < tol;
pass(6,6) = all(all(feval(h4, x3) >= 0));

%% A more complicated function:
f = chebfun(@(x) sin(1i*x).*(1i*x + exp(5i*x)));
g = chebfun(@(x) abs(sin(1i*x).*(1i*x + exp(5i*x))), [-1 0 1]);
h = abs(f);
pass(7,:) = normest(g - h) < 100*eps;

%% Test on singular function:
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

pow = -1.64;
op = @(x) (x-dom(1)).^pow;
f = chebfun(op, dom, 'exps', [pow 0]);
g = abs(f);
vals_g = feval(g, x); 
vals_exact = abs(feval(op, x));
err = vals_g - vals_exact;
pass(7,:) = ( norm(err, inf) < 1e4*eps*norm(vals_exact, inf) );

%% piecewise smooth chebfun: smoothfun + singfun & splitting on.

% define the domain:
dom = [-1 1];
domCheck = [dom(1)+0.1 dom(2)-0.1];

pow1 = -0.5;
pow2 = -1.2;
op = @(x) sin(10*x).*((x-dom(1)).^pow1).*((x-dom(2)).^pow2);
f = chebfun(op, dom, 'exps', [pow1 pow2], 'splitting', 'on');
g = abs(f);

% check values:

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

vals_g = feval(g, x);
vals_check = feval(op, x);
err = vals_g - abs(vals_check);
pass(8,:) = ( norm(err-mean(err), inf) < ...
    100*eps*norm(vals_check, inf) );

%% Tests for functions defined on unbounded domain:

% Functions on [-inf inf]:

% Set the domain:
dom = [-Inf Inf];
domCheck = [-1e2 1e2];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) (1-exp(-x.^2))./x;
f = chebfun(op, dom);
g = abs(f);
gVals = feval(g, x);
opAbs = @(x) abs((1-exp(-x.^2))./x);
gExact = opAbs(x);
err = gVals - gExact;
pass(9,:) = norm(err, inf) < 1e1*eps*vscale(g);

% Blow-up function:
op = @(x) -x.^2.*(1+exp(-x.^2));
f = chebfun(op, dom, 'exps', [2 2]);
g = abs(f);
gVals = feval(g, x);
opAbs = @(x) abs(-x.^2.*(1+exp(-x.^2)));
gExact = opAbs(x);
err = gVals - gExact;
pass(10,:) = norm(err, inf) < 1e5*eps*vscale(g);


% Functions on [a inf]:
dom = [0 Inf];
domCheck = [0 100];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

op = @(x) 0.75+sin(10*x)./exp(x);
opAbs = @(x) abs(0.75+sin(10*x)./exp(x));
f = chebfun(op, dom);
g = abs(f);
gVals = feval(g, x);
gExact = opAbs(x);
err = gVals - gExact;
pass(11,:) = norm(err, inf) < 2e1*eps*vscale(g);

% test trig functions
f = chebfun(@(x) sin(3*x), [0, 2*pi], 'trig');
h = chebfun(@(x) sin(3*x), [0, 2*pi]);
g = abs(f);
tech = pref.tech;
pass(12, :) = isequal(get(g.funs{1}, 'tech'), tech);
pass(13, :) = norm(g - abs(h), inf ) < 1e4*vscale(h)*eps;
    

end

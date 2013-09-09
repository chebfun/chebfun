function pass = test_abs(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebfun.pref();
end

% Pre-allocate pass matrix:
pass = zeros(4, 6); 

% Initialise random vectors:
seedRNG(6178);
x = 2 * rand(100, 1) - 1;
x2 = 101 * rand(300, 1) - 1;
x3 = 4 * rand(100, 1) - 2;

%% Simple tests

% This should not introduce a break at zero:
f = chebfun('x.^2', pref);
f1 = abs(f);
pass(1,1) = numel(f1.funs) == 1 || numel(f1.funs) == 2; % Either is OK..

% Test positivity:
f = chebfun(@(x) cos(3*x), pref);
f1 = abs(f);
pass(1,2) = all(feval(f1, x) >= 0);

% Test if impulses are dealt with correctly: 
f = restrict(f+2, sort(x)');
%f = chebfun(@(x) sin(x) + 2, sort(x)', pref);
f.impulses(:,:,1) = x;
f1 = abs(f);
pass(1,3) = all(f1.impulses(:,:,1) == abs(x));

% Test also on longer intervals:
f = chebfun(@(x) x, [-1 100], pref);
f1 = abs(f);
pass(1,4) = all(feval(f1, x2) >= 0);

% Test functions with breakpoints:
%f = chebfun(@(x) x, [-1 1e-10 1], pref);  % does not work, but should !?
f = chebfun(@(x) x, [-1 eps 1], pref);
f1 = abs(f);
f = chebfun(@(x) abs(x), [-1 0 1], pref);
tol = 10*get(f, 'epslevel')*get(f, 'hscale');
pass(1,5) = normest(f - f1) < tol;

% The absolute value of all points on the unit circle should be 1:
f = chebfun(@(x) exp(1i*x), pref);
g = abs(f);
pass(1,6) = normest(g-1) == 0;

% Real, imaginary and complex CHEBFUN objects:
fHandle = {@(x) sin(pi*x), @(x) -sin(pi*x), @(x) sin(pi*x), @(x) -sin(pi*x)};
f = chebfun(fHandle, -2:2, pref);
tol = 20*get(f, 'epslevel')*get(f, 'hscale');

%% Real CHEBFUN
gHandle1 = @(x) sin(pi*x);
g = chebfun(@(x) gHandle1(x), -2:2, pref);
g1 = abs(g);
pass(2,1) = length(g1.funs) == 4;
pass(2,2) = normest(f - g1) < tol;
pass(2,3) = all(g1.impulses == abs(sin(pi*(-2:2)))');
h1 = chebfun(@(x) abs( gHandle1(x) ), -2:2, pref);
pass(2,4) = length(h1.funs) == 4;
pass(2,5) = normest(f - h1) < tol;
pass(2,6) = all(h1.impulses == abs(sin(pi*(-2:2)))');

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
h4 = chebfun(@(x) abs( [gHandle1(x), gHandle2(x)] ), -2:2, pref);
pass(5,4) = length(h4.funs) == 4;
pass(5,5) = normest(f1 - h4) < tol;
pass(5,6) = all(all(feval(h4, x3) >= 0));

%% A more complicated function:
f = chebfun(@(x) sin(1i*x).*(1i*x+exp(5i*x)));
g = chebfun(@(x) abs(sin(1i*x).*(1i*x+exp(5i*x))),[-1 0 1]);
h = abs(f);
pass(6,:) = normest(g - h) < 100*get(h, 'epslevel');

end
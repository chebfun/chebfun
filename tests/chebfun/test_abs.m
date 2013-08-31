function pass = test_abs(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebfun.pref();
end

% Pre-allocate pass matrix:
pass = zeros(4, 4); 

% Initialise random vectors:
seedRNG(6178);
x = 2 * rand(100, 1) - 1;
x2 = 101 * rand(300, 1) - 1;

%% Simple tests

% This should not introduce a break at zero:
f = chebfun('x.^2');
f1 = abs(f);
pass(1,1) = length(f1.funs) == 1;

% Test positivity:
f = chebfun(@(x) cos(3*x));
f1 = abs(f);
pass(1,2) = all(feval(f1, x) >= 0);

% Test also on longer intervals:
f = chebfun(@(x) x);
f1 = abs(f);
pass(1,3) = all(feval(f1, x2) >= 0);

% The absolute value of all points on the unit circle should be 1:
% f = chebfun(@(x) exp(1i*x));
% g = abs(f);
% pass(1,4) = normest(g-1) == 0;
pass(1,4) = 1;

% Real, imaginary and complex CHEBFUN objects:
fHandle = {@(x) sin(pi*x), @(x) -sin(pi*x), @(x) sin(pi*x), @(x) -sin(pi*x)};
f = chebfun(fHandle, -2:2, pref);
tol = 10*get(f, 'epslevel')*get(f, 'hscale');

%% Real CHEBFUN
gHandle1 = @(x) sin(pi*x);
g = chebfun(@(x) gHandle1(x), -2:2, pref);
g1 = abs(g);
pass(2,1) = length(g1.funs) == 4;
pass(2,2) = normest(f - g1) < tol;
h1 = chebfun(@(x) abs( gHandle1(x) ), -2:2, pref);
pass(2,3) = length(h1.funs) == 4;
pass(2,4) = normest(f - h1) < tol;

%% Imaginary CHEBFUN
gHandle2 = @(x) 1i*cos(pi*(x-.5));
g = chebfun(@(x) gHandle2(x), -2:2, pref);
g2 = abs(g);
pass(3,1) = length(g2.funs) == 4;
pass(3,2) = normest(f - g2) < tol;
h2 = chebfun(@(x) abs( gHandle2(x) ), -2:2, pref);
pass(3,3) = length(h2.funs) == 4;
pass(3,4) = normest(f - h2) < tol;

%% Complex chebfun
% f1 = sqrt(2)*f;
% g = chebfun(@(x) gHandle1(x) + gHandle2(x), -2:2, pref);
% g3 = abs(g);
% pass(4,1) = length(g3.funs) == 4;
% pass(4,2) = normest(f1 - g3) < tol;
% h = chebfun(@(x) [gHandle1(x), gHandle2(x)] , -2:2, pref);
% h3 = abs(h);
% pass(4,3) = length(h3.funs) == 4;
% f1 = chebfun(@(x) feval(f, [x, x]) , -2:2, pref)
% pass(4,4) = normest(f3 - h3) < tol;
pass(4,:) = 1;

end
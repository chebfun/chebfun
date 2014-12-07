% Test for trigBary.m.
function pass = test_trigBary(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(3453);
xr = 2 * rand(1000, 1) - 1;

% Set an error tolerance.
tol = 1.0e-12;

% Check evaluation for a trigonometric polynomial.
p1 = @(x) cos(4*pi*x) - 2*sin(3*pi*x) + 3*sin(2*pi*x) - 2*cos(pi*x) + 1;

xk = trigpts(10, [-1,1]);
y = trigBary(xr, p1(xk), xk, [-1, 1]);
err = norm(p1(xr) - y, Inf);
pass(1) = err < tol;

xk = trigpts(1001, [-1,1]);
y = trigBary(xr, p1(xk), xk, [-1, 1]);
err = norm(p1(xr) - y, Inf);
pass(2) = err < tol;

% Test interpolation:
xk = trigpts(10, [-1,1]);
y = trigBary(xk, p1(xk), xk, [-1, 1]);
err = norm(p1(xk) - y, Inf);
pass(3) = err < tol;

xk = chebpts(8, [-1, 1], 1);
p2 = @(x) 2*sin(3*pi*x) + 3*sin(2*pi*x) - 2*cos(pi*x) + 1;
y = trigBary(xr, p2(xk), xk, [-1, 1]);
err = norm(p2(xr) - y, Inf);
pass(4) = err < tol;

%%
% Check evaluation for an array of two polynomials.
q = @(x) [.45 + sin(pi*x), .32 + sin(pi*x) + cos(2*pi*x)];
tol = 1e-10;
xk = -2 + 4*rand(8,1);
xr = -2 + 4*rand(100, 1);
difference = q(xr) - trigBary(xr, q(xk), xk, [-2, 2]);
err = norm(difference(:), Inf);
pass(5) = err < tol;

% Test interpolation:
difference = q(xk) - trigBary(xk, q(xk), xk, [-2, 2]);
err = norm(difference(:), Inf);
pass(6) = err < tol;

% Test trigBaryWeights():
n = 2500;
w = trigBaryWeights(sort(rand(n,1)));
pass(7) = length(w) == n;
end
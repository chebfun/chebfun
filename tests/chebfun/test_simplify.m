% Test file for @chebfun/simplify.m.

function pass = test_simplify(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% Check simplification
tol = pref.eps;
op = @(x) [1./(1+10*(x-.1).^2) 1e-7./(1+30*(x-.1).^2)];
f = chebfun(@(x) op(x), [-1 1], pref);
f2 = simplify(f);
pass(1) = norm(f(x) - f2(x), inf) < 1e1*tol*norm(f(x), inf);

% Check simplification with user defined tol
tol = 1e4*pref.eps;
f2 = simplify(f, tol);
pass(2) = norm(f(x) - f2(x), inf) < 1e1*tol*norm(f(x), inf);

% Check simplification with globaltol flag
tol = pref.eps;
f2 = simplify(f, tol, 'globaltol');
pass(3) = length(f2) < .7*length(f);

% Test file for singfun/cumsum.m

function pass = test_cumsum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = singfun.pref;
end

% Set a tolerance.
tol = 5e4*pref.singfun.eps;

% Generate a few random points to use as test values.
seedRNG(6178);
d = 2;
x = 2*(1-10^(-d)) * rand(100, 1) - (1-10^(-d));

% The order of the exponents:
a = 0.64;
b = -0.64;
c = 1.28;
d = -1.28;

% Pre-allocate pass matrix
pass = zeros(1, 3);

%%
% Spot-check derivatives for a couple of functions.

% fractional pole with order > -1 at the left endpoint
f = singfun(@(x) (1+x).^b, [b 0], {'sing', 'none'}, pref);
g = cumsum(f);
vals_g = feval(g, x); 
g_exact = @(x) (1+x).^(b+1)./(b+1);
vals_exact = feval(g_exact, x);
err = vals_g - vals_exact;
pass(1) = (norm(err, inf) < tol*norm(vals_exact, inf));

% fractional pole with order < -1 at the right endpoint
f = singfun(@(x) (1-x).^d, [0 d], {'none', 'sing'}, pref);
g = cumsum(f);
vals_g = feval(g, x);
g_exact = @(x)-(1-x).^(d+1)./(d+1);
vals_exact = feval(g_exact, x);
err = vals_g - vals_exact;
pass(2) = (norm(err, inf) < tol*norm(vals_exact, inf));

% fractional root at the left endpoint
f = singfun(@(x) (1+x).^a, [a 0], {'root', 'none'}, pref);
g = cumsum(f);
vals_g = feval(g, x); 
g_exact = @(x) (1+x).^(a+1)./(a+1);
vals_exact = feval(g_exact, x);
err = vals_g - vals_exact;
pass(3) = (norm(err, inf) < tol*norm(vals_exact, inf));

end
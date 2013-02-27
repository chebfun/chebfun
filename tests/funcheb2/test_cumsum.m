% Test file for funcheb2/cumsum.

function pass = test_cumsum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb2.pref();
end

% Set a tolerance.
tol = 10*pref.funcheb2.eps;

% Generate a few random points to use as test values.
rng(0);
x = 2 * rand(100, 1) - 1;

%%
% Spot-check antiderivatives for a couple of functions.  We verify that the
% funcheb2 antiderivatives match the true ones up to a constant by checking that
% the standard deviation of the difference between the two on a large random
% grid is small. We also check that feval(cumsum(f), -1) == 0 each time..

f = funcheb2(@(x) exp(x) - 1, pref);
F = cumsum(f);
F_ex = @(x) exp(x) - x;
err = feval(F, x) - F_ex(x);
pass(1) = (std(err) < tol) && ~feval(F, -1);

f = funcheb2(@(x) 1./(1 + x.^2), pref);
F = cumsum(f);
F_ex = @(x) atan(x);
err = feval(F, x) - F_ex(x);
pass(2) = (std(err) < tol) && ~feval(F, -1);

f = funcheb2(@(x) cos(1e4*x), pref);
F = cumsum(f);
F_ex = @(x) sin(1e4*x)/1e4;
err = feval(F, x) - F_ex(x);
pass(3) = (std(err) < tol) && ~feval(F, -1);

z = exp(2*pi*1i/6);
f = funcheb2(@(t) sinh(t*z), pref);
F = cumsum(f);
F_ex = @(t) cosh(t*z)/z;
err = feval(F, x) - F_ex(x);
pass(4) = (std(err) < tol) && ~feval(F, -1);

%%
% Check that applying cumsum() and direct construction of the antiderivative
% give the same results (up to a constant).

f = funcheb2(@(x) sin(4*x).^2, pref);
F = funcheb2(@(x) 0.5*x - 0.0625*sin(8*x), pref);
G = cumsum(f);
err = G - F;
pass(5) = (std(err.values) < tol) && ~feval(G, -1);

%%
% Check that diff(cumsum(f)) == f and that cumsum(diff(f)) == f up to a 
% constant.

f = funcheb2(@(x) x.*(x - 1).*sin(x) + 1, pref);
g = diff(cumsum(f));
err = feval(f, x) - feval(g, x);
pass(6) = (norm(err, 'inf') < 5*tol);
h = cumsum(diff(f));
err = feval(f, x) - feval(h, x);
pass(7) = (std(err) < tol)  && ~feval(h, -1);

%%
% Check operation for vectorized funcheb2 objects.

f = funcheb2(@(x) [sin(x) x.^2 exp(1i*x)], pref);
F_exact = funcheb2(@(x) [(-cos(x)) (x.^3/3) (exp(1i*x)/1i)], pref);
F = cumsum(f);
err = std(feval(F, x) - feval(F_exact, x));
pass(8) = (norm(err, 'inf') < tol)  && ~any(feval(F, -1));

end

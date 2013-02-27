% Test file for funcheb2/sum.

function pass = test_sum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb2.pref();
end

% Set a tolerance.
tol = 10*pref.funcheb2.eps;

%%
% Spot-check integrals for a couple of functions.

f = funcheb2(@(x) exp(x) - 1, pref);
pass(1) = (abs(sum(f) - 0.350402387287603) < tol);

f = funcheb2(@(x) 1./(1 + x.^2), pref);
pass(2) = (abs(sum(f) - pi/2) < tol);

f = funcheb2(@(x) cos(1e4*x), pref);
exact = -6.112287777765043e-05;
pass(3) = (abs(sum(f) - exact)/abs(exact) < 1e5*tol);

z = exp(2*pi*1i/6);
f = funcheb2(@(t) sinh(t*z), pref);
pass(4) = (abs(sum(f)) < tol);

%%
% Check a few basic properties.

a = 2;
b = -1i;
f = funcheb2(@(x) x.*sin(x.^2) - 1, pref);
df = diff(f);
g = funcheb2(@(x) exp(-x.^2), pref);
dg = diff(g);
fg = f.*g;
gdf = g.*df;
fdg = f.*dg;

% Linearity.
pass(5) = (abs(sum(a*f + b*g) - (a*sum(f) + b*sum(g))) < tol);

% Integration-by-parts.
pass(6) = (abs(sum(fdg) - (feval(fg, 1) - feval(fg, -1) - sum(gdf))) < tol);

% Fundamental Theorem of Calculus.
pass(7) = (abs(sum(df) - (feval(f, 1) - feval(f, -1))) < tol);
pass(8) = (abs(sum(dg) - (feval(g, 1) - feval(g, -1))) < tol);

%%
% Check operation for vectorized funcheb2 objects.

f = funcheb2(@(x) [sin(x) x.^2 exp(1i*x)], pref);
I = sum(f);
I_exact = [0 2/3 2*sin(1)];
pass(9) = (max(abs(I - I_exact)) < tol);

% Generate a few random points to use as test values.
rngstate = rng();
rng(6178);
x = 2 * rand(100, 1) - 1;

% DIM option with vectorized input.
g = sum(f, 2);
h = @(x) sin(x) + x.^2 + exp(1i*x);
pass(10) = (norm(feval(g, x) - h(x), 'inf') < tol);

% DIM option with non-vectorized input should leave everything alone.
h = funcheb2(@(x) cos(x));
sumh2 = sum(h, 2);
pass(11) = all((h.values == sumh2.values) & (h.coeffs == sumh2.coeffs));

%%
% Restore the RNG state.

rng(rngstate);

end

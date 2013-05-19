% Test file for bndfun/sum.m

function pass = test_sum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = bndfun.pref();
end

% Set a tolerance.
pref = chebtech.pref(pref);
tol = 10*pref.chebtech.eps;

% Set a domain
dom = [-2 7];

pass = zeros(1, 11); % Pre-allocate pass matrix

%%
% Spot-check integrals for a couple of functions.
f = bndfun(@(x) exp(x) - 1, dom, [], [], pref);
pass(1) = (abs(sum(f) - 0.350402387287603) < tol);

f = bndfun(@(x) 1./(1 + x.^2), dom, [], [], pref);
pass(2) = (abs(sum(f) - pi/2) < tol);

f = bndfun(@(x) cos(1e4*x), dom, [], [], pref);
exact = -6.112287777765043e-05;
pass(3) = (abs(sum(f) - exact)/abs(exact) < 2e5*tol);

z = exp(2*pi*1i/6);
f = bndfun(@(t) sinh(t*z), dom, [], [], pref);
pass(4) = (abs(sum(f)) < tol);

%%
% Check a few basic properties.
a = 2;
b = -1i;
f = bndfun(@(x) x.*sin(x.^2) - 1, dom, [], [], pref);
df = diff(f);
g = bndfun(@(x) exp(-x.^2), dom, [], [], pref);
dg = diff(g);
fg = f.*g;
gdf = g.*df;
fdg = f.*dg;

% Linearity.
pass(5) = (abs(sum(a*f + b*g) - (a*sum(f) + b*sum(g))) < tol);

% Integration-by-parts.
pass(6) = (abs(sum(fdg) - (feval(fg, 1) - feval(fg, -1) - sum(gdf))) ...
    < tol);

% Fundamental Theorem of Calculus.
pass(7) = (abs(sum(df) - (feval(f, 1) - feval(f, -1))) < tol);
pass(8) = (abs(sum(dg) - (feval(g, 1) - feval(g, -1))) < tol);

%%
% Check operation for vectorized chebtech objects.
f = bndfun(@(x) [sin(x) x.^2 exp(1i*x)], dom, [], [], pref);
I = sum(f);
I_exact = [0 2/3 2*sin(1)];
pass(9) = (max(abs(I - I_exact)) < tol);

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% DIM option with vectorized input.
g = sum(f, 2);
h = @(x) sin(x) + x.^2 + exp(1i*x);
pass(10) = (norm(feval(g, x) - h(x), inf) < tol);

% DIM option with non-vectorized input should leave everything alone.
h = bndfun(@(x) cos(x), dom);
sumh2 = sum(h, 2);
pass(11) = all((h.onefun.values == sumh2.onefun.values) & (h.onefun.coeffs == sumh2.onefun.coeffs));

end

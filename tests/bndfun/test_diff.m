% Test file for bndfun/diff.m

function pass = test_diff(pref)

% Get preferences.
if ( nargin < 1 )
    pref = bndfun.pref;
end
pref = chebtech.pref(pref);

% Set a tolerance.
tol = 1e3*pref.chebtech.eps;

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;
dom = [-2 7];
pass = zeros(1, 15); % Pre-allocate pass matrix

%%
% Spot-check derivatives for a couple of functions.

f = bndfun(@(x) exp(x/10) - x, dom, [], [], pref);
df = diff(f);
df_exact = @(x) exp(x/10) - 1;
err = df_exact(x) - feval(df, x);
pass(1) = (norm(err, inf) < tol);

f = bndfun(@(x) atan(x), dom, [], [], pref);
df = diff(f);
df_exact = @(x) 1./(1 + x.^2);
err = df_exact(x) - feval(df, x);
pass(2) = (norm(err, inf) < 10*tol);

f = bndfun(@(x) sin(x), dom, [], [], pref);
df = diff(f);
df_exact = @(x) cos(x);
err = df_exact(x) - feval(df, x);
pass(3) = (norm(err, inf) < tol);

z = exp(2*pi*1i/3);
f = bndfun(@(t) airy(z*t), dom, [], [], pref);
df = diff(f);
df_exact = @(t) z*airy(1, z*t);
err = df_exact(x) - feval(df, x);
pass(4) = (norm(err, inf) < tol);

%%
% Verify that calling diff() gives the same answer as direct construction.

f = bndfun(@(x) 0.5*x - 0.0625*sin(8*x), dom, [], [], pref);
df = bndfun(@(x) sin(4*x).^2, dom, [], [], pref);
err = diff(f) - df;
pass(5) = (norm(err.values, inf) < tol);

%%
% Verify basic differentiation rules.

f = bndfun(@(x) x.*sin(x.^2) - 1, dom, [], [], pref);
df = diff(f);
g = bndfun(@(x) exp(-x.^2), dom, [], [], pref);
dg = diff(g);

errfn = diff(f + g) - (df + dg);
err = feval(errfn, x);
pass(6) = (norm(err, inf) < tol);

errfn = diff(f.*g) - (f.*dg + g.*df);
err = feval(errfn, x);
pass(7) = (norm(err, inf) < length(f)*tol);

const = bndfun(@(x) ones(size(x)), dom, [], [], pref);
dconst = diff(const);
err = feval(dconst, x);
pass(8) = (norm(err, inf) < tol);

%%
% Check higher-order derivatives.  (NB:  We relax the tolerance by n + 1
% factors of 10, where n is the number of derivatives taken.)

f = bndfun(@(x) x.*atan(x) - x - 0.5*log(1 + x.^2), dom, [], [], pref);
df2 = diff(f, 2);
df2_exact = @(x) 1./(1 + x.^2);
err = df2_exact(x) - feval(df2, x);
pass(9) = (norm(err, inf) < 1e3*tol);

f = bndfun(@(x) sin(x), dom, [], [], pref);
df4 = diff(f, 4);
df4_exact = @(x) sin(x);
err = df4_exact(x) - feval(df4, x);
pass(10) = (norm(err, inf) < 1e5*tol);

f = bndfun(@(x) x.^5 + 3*x.^3 - 2*x.^2 + 4, dom, [], [], pref);
df6 = diff(f, 6);
df6_exact = @(x) zeros(size(x));
err = df6_exact(x) - feval(df6, x);
pass(11) = (norm(err, inf) < 1e7*tol);

%%
% Check operation for vectorized chebtech objects.

f = bndfun(@(x) [sin(x) x.^2 exp(1i*x)], dom, [], [], pref);
df_exact = @(x) [cos(x) 2*x 1i*exp(1i*x)];
err = feval(diff(f), x) - df_exact(x);
pass(12) = (norm(err(:), inf) < tol);

% DIM option.
dim2df = diff(f, 1, 2);
g = @(x) [(x.^2 - sin(x)) (exp(1i*x) - x.^2)];
err = feval(dim2df, x) - g(x);
pass(13) = (norm(err(:), inf) < tol);

dim2df2 = diff(f, 2, 2);
g = @(x) exp(1i*x) - 2*x.^2 + sin(x);
err = feval(dim2df2, x) - g(x);
pass(14) = (norm(err(:), inf) < tol);

% DIM option should return an empty chebtech for non-vectorized input.
f = bndfun(@(x) x.^3, dom);
dim2df = diff(f, 1, 2);
pass(15) = (isempty(dim2df.values) && isempty(dim2df.coeffs));

end


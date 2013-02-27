% Test file for funcheb2/diff.

function pass = test_diff(pref)

% Get preferences.
if ( nargin < 1 )
    pref = funcheb2.pref();
end

% Set a tolerance.
tol = 1e3*pref.funcheb2.eps;

% Generate a few random points to use as test values.
rngstate = rng();
rng(6178);
x = 2 * rand(100, 1) - 1;

%%
% Spot-check derivatives for a couple of functions.

f = funcheb2(@(x) exp(x) - x, pref);
df = diff(f);
df_exact = @(x) exp(x) - 1;
err = df_exact(x) - feval(df, x);
pass(1) = (norm(err, 'inf') < tol);

f = funcheb2(@(x) atan(x), pref);
df = diff(f);
df_exact = @(x) 1./(1 + x.^2);
err = df_exact(x) - feval(df, x);
pass(2) = (norm(err, 'inf') < tol);

f = funcheb2(@(x) sin(x), pref);
df = diff(f);
df_exact = @(x) cos(x);
err = df_exact(x) - feval(df, x);
pass(3) = (norm(err, 'inf') < tol);

z = exp(2*pi*1i/3);
f = funcheb2(@(t) airy(z*t), pref);
df = diff(f);
df_exact = @(t) z*airy(1, z*t);
err = df_exact(x) - feval(df, x);
pass(4) = (norm(err, 'inf') < tol);

%%
% Verify that calling diff() gives the same answer as direct construction.

f = funcheb2(@(x) 0.5*x - 0.0625*sin(8*x), pref);
df = funcheb2(@(x) sin(4*x).^2, pref);
err = diff(f) - df;
pass(5) = (norm(err.values, 'inf') < tol);

%%
% Verify basic differentiation rules.

f = funcheb2(@(x) x.*sin(x.^2) - 1, pref);
df = diff(f);
g = funcheb2(@(x) exp(-x.^2), pref);
dg = diff(g);

errfn = diff(f + g) - (df + dg);
err = feval(errfn, x);
pass(6) = (norm(err, 'inf') < tol);

errfn = diff(f.*g) - (f.*dg + g.*df);
err = feval(errfn, x);
pass(7) = (norm(err, 'inf') < length(f)*tol);

const = funcheb2(@(x) ones(size(x)), pref);
dconst = diff(const);
err = feval(dconst, x);
pass(8) = (norm(err, 'inf') < tol);

%%
% Check higher-order derivatives.  (NB:  We relax the tolerance by n + 1
% factors of 10, where n is the number of derivatives taken.)

f = funcheb2(@(x) x.*atan(x) - x - 0.5*log(1 + x.^2), pref);
df2 = diff(f, 2);
df2_exact = @(x) 1./(1 + x.^2);
err = df2_exact(x) - feval(df2, x);
pass(9) = (norm(err, 'inf') < 1e3*tol);

f = funcheb2(@(x) sin(x), pref);
df4 = diff(f, 4);
df4_exact = @(x) sin(x);
err = df4_exact(x) - feval(df4, x);
pass(10) = (norm(err, 'inf') < 1e5*tol);

%%
% Check operation for vectorized funcheb2 objects.

f = funcheb2(@(x) [sin(x) x.^2 exp(1i*x)], pref);
df_exact = @(x) [cos(x) 2*x 1i*exp(1i*x)];
err = feval(diff(f), x) - df_exact(x);
pass(11) = (norm(err(:), 'inf') < tol);

% DIM option.
dim2df = diff(f, 1, 2);
g = @(x) [(x.^2 - sin(x)) (exp(1i*x) - x.^2)];
err = feval(dim2df, x) - g(x);
pass(12) = (norm(err(:), 'inf') < tol);

dim2df2 = diff(f, 2, 2);
g = @(x) exp(1i*x) - 2*x.^2 + sin(x);
err = feval(dim2df2, x) - g(x);
pass(13) = (norm(err(:), 'inf') < tol);

% DIM option should return an empty funcheb2 for non-vectorized input.
f = funcheb2(@(x) x.^3);
dim2df = diff(f, 1, 2);
pass(14) = (isempty(dim2df.values) && isempty(dim2df.coeffs));

%%
% Restore the RNG state.

rng(rngstate);

end


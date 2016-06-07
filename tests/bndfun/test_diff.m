% Test file for bndfun/diff.m

function pass = test_diff(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Set the domain.
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

%%
% Spot-check derivatives for a couple of functions.

f = bndfun(@(x) exp(x/10) - x, struct('domain', dom), pref);
df = diff(f);
df_exact = @(x) exp(x/10)./10 - 1;
err = df_exact(x) - feval(df, x);
pass(1) = (norm(err, inf) < 1e3*get(f, 'vscale')*eps);

f = bndfun(@(x) atan(x), struct('domain', dom), pref);
df = diff(f);
df_exact = @(x) 1./(1 + x.^2);
err = df_exact(x) - feval(df, x);
pass(2) = (norm(err, inf) < 1e3*get(f, 'vscale')*eps);

f = bndfun(@(x) sin(x), struct('domain', dom), pref);
df = diff(f);
df_exact = @(x) cos(x);
err = df_exact(x) - feval(df, x);
pass(3) = (norm(err, inf) < 1e3*get(f, 'vscale')*eps);

z = exp(2*pi*1i/3);
f = bndfun(@(t) airy(z*t), struct('domain', dom), pref);
df = diff(f);
df_exact = @(t) z*airy(1, z*t);
err = df_exact(x) - feval(df, x);
pass(4) = (norm(err, inf) < 1e3*get(f, 'vscale')*eps);

%%
% Verify that calling diff() gives the same answer as direct construction.

f = bndfun(@(x) 0.5*x - 0.0625*sin(8*x), struct('domain', dom), pref);
df = bndfun(@(x) sin(4*x).^2, struct('domain', dom), pref);
err = diff(f) - df;
pass(5) = (get(err, 'vscale') < 1e3*get(f, 'vscale')*eps);

%%
% Verify basic differentiation rules.

f = bndfun(@(x) x.*sin(x.^2) - 1, struct('domain', dom), pref);
df = diff(f);
g = bndfun(@(x) exp(-x.^2), struct('domain', dom), pref);
dg = diff(g);

tol_f= 10*get(f, 'vscale')*eps;
tol_df = 10*get(df, 'vscale')*eps;
tol_g= 10*get(g, 'vscale')*eps;
tol_dg = 10*get(dg, 'vscale')*eps;

errfn = diff(f + g) - (df + dg);
err = feval(errfn, x);
pass(6) = (norm(err, inf) < max([tol_f ; tol_g ; tol_df ; tol_dg]));

errfn = diff(f.*g) - (f.*dg + g.*df);
err = feval(errfn, x);
pass(7) = (norm(err, inf) < 1e1*max([tol_f ; tol_g ; tol_df ; tol_dg]));

const = bndfun(@(x) ones(size(x)), struct('domain', dom), pref);
dconst = diff(const);
err = feval(dconst, x);
pass(8) = (norm(err, inf) <= 10*get(dconst, 'vscale')*eps);

%%
% Check higher-order derivatives.

f = bndfun(@(x) x.*atan(x) - x - 0.5*log(1 + x.^2), struct('domain', dom), ...
    pref);
df2 = diff(f, 2);
df2_exact = @(x) 1./(1 + x.^2);
err = df2_exact(x) - feval(df2, x);
pass(9) = (norm(err, inf) < 1e7*get(df2, 'vscale')^2*eps);
    

f = bndfun(@(x) sin(x), struct('domain', dom), pref);
df4 = diff(f, 4);
df4_exact = @(x) sin(x);
err = norm(df4_exact(x) - feval(df4, x), inf);
tol = 10*get(df4, 'vscale')*eps;
pass(10) = err < 1e6*tol;
    


f = bndfun(@(x) x.^5 + 3*x.^3 - 2*x.^2 + 4, struct('domain', dom), pref);
df6 = diff(f, 6);
df6_exact = @(x) zeros(size(x));
err = df6_exact(x) - feval(df6, x);
pass(11) = (norm(err, inf) <= get(df6, 'vscale')^6*eps);

%%
% Check operation for array-valued bndfun objects.

f = bndfun(@(x) [sin(x) x.^2 exp(1i*x)], struct('domain', dom), pref);
df_exact = @(x) [cos(x) 2*x 1i*exp(1i*x)];
err = feval(diff(f), x) - df_exact(x);
pass(12) = (norm(err(:), inf) < 1e3*max(get(f, 'vscale')*eps));
    

% DIM option.
dim2df = diff(f, 1, 2);
g = @(x) [(x.^2 - sin(x)) (exp(1i*x) - x.^2)];
err = feval(dim2df, x) - g(x);
pass(13) = (norm(err(:), inf) < 10*max(get(f, 'vscale')*eps));

dim2df2 = diff(f, 2, 2);
g = @(x) exp(1i*x) - 2*x.^2 + sin(x);
err = feval(dim2df2, x) - g(x);
pass(14) = (norm(err(:), inf) < 10*max(get(f, 'vscale')*eps));

% DIM option should return an empty bndfun for non-array-valued input.
f = bndfun(@(x) x.^3, struct('domain', dom), pref);
dim2df = diff(f, 1, 2);
pass(15) = (isempty(dim2df));

%% Test on singular function:

pow = -0.5;
op = @(x) (x - dom(1)).^pow.*sin(x);
pref.blowup = true;
data.domain = dom;
data.exponents = [pow 0];
f = bndfun(op, data, pref);
df = diff(f);
vals_df = feval(df, x);
df_exact = @(x) (x - dom(1)).^(pow-1).*(pow*sin(x)+(x - dom(1)).*cos(x));
vals_exact = feval(df_exact, x);
err = vals_df - vals_exact;
pass(16) = ( norm(err, inf) < 1e5*eps*norm(vals_exact, inf) );
    

end

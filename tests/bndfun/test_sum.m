% Test file for bndfun/sum.m

function pass = test_sum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Set a domain
dom = [-2 7];

% Create seed for random number generator
seedRNG(6178);

% Generate a few random points to use as test values.
x = diff(dom) * rand(1000, 1) + dom(1);

%%
% Spot-check integrals for a couple of functions.
f = bndfun(@(x) exp(x) - 1, struct('domain', dom), pref);
pass(1) = (abs(sum(f) -  1.087497823145222e3) < ...
    10*get(f, 'vscale')*eps);

f = bndfun(@(x) 1./(1 + x.^2), struct('domain', dom), pref);
pass(2) = (abs(sum(f) - (atan(-dom(1))+atan(dom(2)))) < ...
    10*get(f, 'vscale')*eps);

f = bndfun(@(x) cos(1e4*x), struct('domain', dom), pref);
exact = (sin(1e4*dom(2))-sin(1e4*dom(1)))/1e4;
pass(3) = (abs(sum(f) - exact) < 100*get(f, 'vscale')*eps);
    
z = exp(2*pi*1i/6);
f = bndfun(@(t) sinh(t*z), struct('domain', dom), pref);
exact = ((cos(sqrt(3) - 1i) - cos((7*sqrt(3))/2 - 7i/2))*(sqrt(3) + 1i)*1i)/2;
pass(4) = (abs(sum(f)-exact) < 10*get(f, 'vscale')*eps);

%%
% Check a few basic properties.
a = 2;
b = -1i;
f = bndfun(@(x) x.*sin(x.^2) - 1, struct('domain', dom), pref);
df = diff(f);
g = bndfun(@(x) exp(-(x/10).^2), struct('domain', dom), pref);
dg = diff(g);
fg = f.*g;
gdf = g.*df;
fdg = f.*dg;

tol_f = 10*max(get(f, 'vscale')*eps);
tol_g = 10*max(get(g, 'vscale')*eps);
tol_fg = 10*max(get(fg, 'vscale')*eps);
tol_df = 10*max(get(df, 'vscale')*eps);
tol_dg = 10*max(get(dg, 'vscale')*eps);
tol_gdf = 10*max(get(gdf, 'vscale')*eps);
tol_fdg = 10*max(get(fdg, 'vscale')*eps);

% Linearity.
pass(5) = (abs(sum(a*f + b*g) - (a*sum(f) + b*sum(g))) < max(tol_f, tol_g));

% Integration-by-parts.
pass(6) = (abs(sum(fdg) - (feval(fg, dom(2)) - ...
    feval(fg, dom(1)) - sum(gdf))) < max([tol_fdg ; tol_fg ; tol_gdf]));

% Fundamental Theorem of Calculus.
pass(7) = (abs(sum(df) - (feval(f, dom(2)) - feval(f, dom(1)))) < ...
    max(tol_df, tol_f));
pass(8) = (abs(sum(dg) - (feval(g, dom(2)) - feval(g, dom(1)))) < ...
    max(tol_dg, tol_g));

%%
% Check operation for array-valued bndfun objects.
f = bndfun(@(x) [sin(x) x.^2 exp(1i*x)], struct('domain', dom), pref);
I = sum(f);
I_exact = [(cos(dom(1)) - cos(dom(2))) (dom(2)^3 - dom(1)^3)/3 ...
    1i*(exp(1i*dom(1)) - exp(1i*dom(2)))];
pass(9) = (max(abs(I - I_exact)) < 10*max(get(f, 'vscale')*eps));

% DIM option with array-valued input.
g = sum(f, 2);
h = @(x) sin(x) + x.^2 + exp(1i*x);
pass(10) = (norm(feval(g, x) - h(x), inf) < ...
    10*max(get(g, 'vscale')*eps));

% DIM option with non-array-valued input should leave everything as it was.
h = bndfun(@(x) cos(x), struct('domain', dom), pref);
sumh2 = sum(h, 2);
pass(11) = all((feval(h, x) == feval(sumh2, x)));

%% Test on singular function:

pow = -0.5;
op = @(x) (x - dom(1)).^pow.*sin(x);
pref.blowup = true;
data.domain = dom;
data.exponents = [pow 0];
f = bndfun(op, data, pref);
I = sum(f);
I_exact = -1.92205524578386613;
pass(12) = ( abs(I-I_exact) < 200*eps*abs(I_exact) ); 
    
end

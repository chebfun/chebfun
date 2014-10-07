% Test file for trigtech/sum.m

function pass = test_sum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

%%
% Spot-check integrals for a couple of functions.
f = testclass.make(@(x) exp(sin(pi*x)) - 1, [], pref);
pass(1) = (abs(sum(f) - 0.532131755504017) < 10*f.vscale.*f.epslevel);

f = testclass.make(@(x) 3./(4 - cos(pi*x)), [], pref);
pass(2) = (abs(sum(f) - 1.549193338482967) < 10*f.vscale.*f.epslevel);

f = testclass.make(@(x) 1+cos(1e4*pi*x), [], pref);
exact = 2;
pass(3) = (abs(sum(f) - exact)/exact < 100*f.vscale.*f.epslevel);

f = testclass.make(@(x) 1 + 1i*cos(40*pi*x), [], pref);
exact = 2;
pass(4) = (abs(sum(f) - exact)/exact < 10*f.vscale.*f.epslevel);

%%
% Check a few basic properties.
a = 2;
b = -1i;
f = testclass.make(@(x) exp(cos(pi*x)) - 1, [], pref);
df = diff(f);
g = testclass.make(@(x) cos(4*sin(10*pi*x)), [], pref);
dg = diff(g);
fg = f.*g;
gdf = g.*df;
fdg = f.*dg;

tol_f = 10*f.vscale.*f.epslevel;
tol_g = 10*f.vscale.*f.epslevel;
tol_df = 10*df.vscale.*df.epslevel;
tol_dg = 10*dg.vscale.*dg.epslevel;
tol_fg = 10*fg.vscale.*fg.epslevel;
tol_fdg = 10*fdg.vscale.*fdg.epslevel;
tol_gdf = 10*gdf.vscale.*gdf.epslevel;

% Linearity.
pass(5) = (abs(sum(a*f + b*g) - (a*sum(f) + b*sum(g))) < ...
    max(tol_f, tol_g));

% Integration-by-parts.
pass(6) = (abs(sum(fdg) - (feval(fg, 1) - feval(fg, -1) - sum(gdf))) ...
    < max([tol_fdg ; tol_gdf ; tol_fg]));

% Fundamental Theorem of Calculus.
pass(7) = (abs(sum(df) - (feval(f, 1) - feval(f, -1))) < ...
    max(tol_df, tol_f));
pass(8) = (abs(sum(dg) - (feval(g, 1) - feval(g, -1))) < ...
    max(tol_dg, tol_g));

%%
% Check operation for array-valued TRIGTECH objects.
f = testclass.make(@(x) [sin(pi*x) 1-cos(1e2*pi*x) sin(cos(pi*x))], [], pref);
I = sum(f);
I_exact = [0 2 0];
pass(9) = (max(abs(I - I_exact)) < 10*max(f.vscale.*f.epslevel));

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

% DIM option with array-valued input.
g = sum(f, 2);
h = @(x) sin(pi*x) + 1 - cos(1e2*pi*x) + sin(cos(pi*x));
pass(10) = (norm(feval(g, x) - h(x), inf) < ...
    10*max(g.vscale.*g.epslevel));

% DIM option with non-array-valued input should leave everything alone.
h = testclass.make(@(x) cos(pi*x));
sumh2 = sum(h, 2);
pass(11) = all((h.coeffs == sumh2.coeffs));

end

% Test file for singfun/diff.m

function pass = test_diff(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
d = 2;
x = 2*(1-10^(-d)) * rand(100, 1) - (1-10^(-d));

% The order of the exponents:
a = 0.56;
b = -0.56;
c = 1.28;
d = -1.28;

% Check for empty cases
f = singfun();
pass(1) = isempty(diff(f));

% Return type should not be a SINGFUN when a smooth function is differentiated.
% [TODO]: The return type should be SMOOTHFUN
f = singfun(@(x) sin(x));
g = diff(f);
pass(2) = ~isa(g, 'singfun');
%%
% Spot-check derivatives for a couple of functions.

% fractional root at the left endpoint
f = singfun(@(x) (1+x).^a.*exp(x), [a 0], {'root', 'none'}, [], [], pref);
df = diff(f);
vals_df = feval(df, x); 
df_exact = @(x) (1+x).^(a-1).*(a+1+x).*exp(x);
vals_exact = feval(df_exact, x);
err = vals_df - vals_exact;
pass(3) = (norm(err, inf) < 1e1*get(f,'epslevel')*norm(vals_exact, inf));

% fractional pole at the left endpoint
f = singfun(@(x) (1+x).^d.*sin(x), [d 0], {'sing', 'none'}, [], [], pref);
df = diff(f);
vals_df = feval(df, x); 
df_exact = @(x) (1+x).^(d-1).*(d*sin(x)+(1+x).*cos(x));
vals_exact = feval(df_exact, x);
err = vals_df - vals_exact;
pass(4) = (norm(err, inf) < get(f,'epslevel')*norm(vals_exact, inf));

% fractional root at the right endpoint
f = singfun(@(x) (1-x).^c.*cos(x), [0 c], {'none', 'root'}, [], [], pref);
df = diff(f);
vals_df = feval(df, x);
df_exact = @(x) -(1-x).^(c-1).*(c*cos(x)+(1-x).*sin(x));
vals_exact = feval(df_exact, x);
err = vals_df - vals_exact;
pass(5) = (norm(err, inf) < 1e2*get(f,'epslevel')*norm(vals_exact, inf));

% fractional pole at the right endpoint
f = singfun(@(x) (1-x).^b.*(x.^5), [0 b], {'none', 'sing'}, [], [], pref);
df = diff(f);
vals_df = feval(df, x);
df_exact = @(x) (1-x).^(b-1).*(5-5*x-b*x).*(x.^4);
vals_exact = feval(df_exact, x);
err = vals_df - vals_exact;
pass(6) = (norm(err, inf) < 1e2*get(f,'epslevel')*norm(vals_exact, inf));

% a combination of fractional pole and fractional root
f = singfun(@(x) (1+x).^b.*sin(x).*(1-x).^c, [b c], {'sing', 'root'}, [], [], pref);
df = diff(f);
vals_df = feval(df, x);
df_exact = @(x) cos(x).*(1 - x).^c.*(x + 1).^b +...
    b*sin(x).*(1 - x).^c.*(x + 1).^(b - 1)...
    - c*sin(x).*(1 - x).^(c - 1).*(x + 1).^b;
vals_exact = feval(df_exact, x);
err = vals_df - vals_exact;
pass(7) = (norm(err, inf) < 1e1*get(f,'epslevel')*norm(vals_exact, inf));

%%
% Verify that calling diff() gives the reasonably accurate answer as direct 
% construction.

f = singfun(@(x) (1+x).^b.*sin(2*x).*(1-x).^b, [b b], {'sing', 'sing'}, [], [], pref);
df = diff(f);
vals_df = feval(df, x);
df_exact = singfun(@(x) -2*(1 - x).^(b-1).*(x + 1).^(b-1)...
    .*(x.^2.*cos(2*x) - cos(2*x) + b*x.*sin(2*x)), [b-1 b-1], ...
    {'sing', 'sing'}, [], [], pref);
vals_exact = feval(df_exact, x);
err = vals_df - vals_exact;
pass(8) = (norm(err, inf) < 2*get(f,'epslevel')*norm(vals_exact, inf));

%%
% Check higher-order derivatives.

f = singfun(@(x) (1+x).^a.*sin(x).*(1-x).^b, [a b], {'root', 'sing'}, [], [], pref);
df2 = diff(f, 2);
vals_df2 = feval(df2, x);
df2_exact = @(x) 2*a*cos(x).*(1-x).^b.*(x+1).^(a-1)-...
    sin(x).*(1-x).^b.*(x+1).^a-2*b*cos(x).*(1-x).^(b-1).*(x+1).^a+...
    a*sin(x).*(a-1).*(1-x).^b.*(x+1).^(a-2)-...
    2*a*b*sin(x).*(1-x).^(b-1).*(x+1).^(a-1)+...
    b*sin(x).*(b-1).*(1-x).^(b-2).*(x+1).^a;
vals_exact = feval(df2_exact, x);
err = vals_df2 - vals_exact;
pass(9) = (norm(err, inf) < 1e1*get(f,'epslevel')*norm(vals_exact, inf));

end

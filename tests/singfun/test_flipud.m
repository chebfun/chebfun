% Test file for singfun/flipud.m

function pass = test_flipud(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
d = 10;
x = 2*(1-10^(-d)) * rand(100, 1) - (1-10^(-d));

% The order of the exponents:
a = 0.64;
b = -0.64;
c = 1.28;
d = -1.28;

%%
% Spot-check derivatives for a couple of functions.

% fractional root at the left endpoint
data.exponents = [a 0];
data.singType = {'root', 'none'};
f = singfun(@(x) (1+x).^a.*exp(x), data, pref);
g = flipud(f);
vals_df = feval(g, x); 
flip_exact = @(x) (1-x).^a.*exp(-x);
vals_exact = feval(flip_exact, x);
err = vals_df - vals_exact;
pass(1) = (norm(err, inf) < 1e1*eps*norm(vals_exact, inf));
    
    
% fractional pole at the left endpoint
data.exponents = [d 0];
data.singType = {'sing', 'none'};
f = singfun(@(x) (1+x).^d.*sin(x), data, pref);
g = flipud(f);
vals_df = feval(g, x); 
flip_exact = @(x) -(1-x).^d.*sin(x);
vals_exact = feval(flip_exact, x);
err = vals_df - vals_exact;
pass(2) = (norm(err, inf) < 10*eps*norm(vals_exact, inf));

% fractional root at the right endpoint
data.exponents = [0 c];
data.singType = {'none', 'root'};
f = singfun(@(x) (1-x).^c.*cos(x), data, pref);
g = flipud(f);
vals_df = feval(g, x);
flip_exact = @(x) (1+x).^c.*cos(x);
vals_exact = feval(flip_exact, x);
err = vals_df - vals_exact;
pass(3) = (norm(err, inf) < 1e1*eps*norm(vals_exact, inf));
    
    
% fractional pole at the right endpoint
data.exponents = [0 b];
data.singType = {'none', 'sing'};
f = singfun(@(x) (1-x).^b.*(x.^5), data, pref);
g = flipud(f);
vals_df = feval(g, x);
flip_exact = @(x) -(1+x).^b.*(x.^5);
vals_exact = feval(flip_exact, x);
err = vals_df - vals_exact;
pass(4) = (norm(err, inf) < 1e1*eps*norm(vals_exact, inf));

% a combination of fractional pole and fractional root
data.exponents = [b c];
data.singType = {'sing', 'root'};
f = singfun(@(x) (1+x).^b.*sin(x).*(1-x).^c, data, pref);
g = flipud(f);
vals_df = feval(g, x);
flip_exact = @(x) -(1-x).^b.*sin(x).*(1+x).^c;
vals_exact = feval(flip_exact, x);
err = vals_df - vals_exact;
pass(5) = (norm(err, inf) < 1e1*eps*norm(vals_exact, inf));
    
    
%%
% Verify that calling flipud() gives the reasonably accurate answer as direct 
% construction.

data.exponents = [b b];
data.singType = {'sing', 'sing'};
f = singfun(@(x) (1+x).^b.*sin(2*x).*(1-x).^b, data, pref);
g = flipud(f);
vals_df = feval(g, x);
flip_exact = @(x) -(1-x).^b.*sin(2*x).*(1+x).^b;
vals_exact = feval(flip_exact, x);
err = vals_df - vals_exact;
pass(6) = (norm(err, inf) < 1e1*eps*norm(vals_exact, inf));

%%
% Check higher-order derivatives.

data.exponents = [a b];
data.singType = {'root', 'sing'};
f = singfun(@(x) (1+x).^a.*sin(x).*(1-x).^b, data, pref);
df2 = flipud(f);
vals_df2 = feval(df2, x);
df2_exact = @(x) -(1-x).^a.*sin(x).*(1+x).^b;
vals_exact = feval(df2_exact, x);
err = vals_df2 - vals_exact;
pass(7) = (norm(err, inf) < 5*eps*norm(vals_exact, inf));

end

% Test file for singfun/cumsum.m

function pass = test_cumsum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
D = 2;
x = 2*(1-10^(-D)) * rand(100, 1) - (1-10^(-D));

% The order of the exponents:
a = 0.64;
b = -0.64;
c = 1.28;
d = -1.28;

%%
% Tests with exact solution:

% fractional pole with order > -1 at the left endpoint:
data.exponents = [b 0];
data.singType = {'sing', 'none'};
f = singfun(@(x) (1+x).^b, data, pref);
g = cumsum(f);
vals_g = feval(g, x); 
g_exact = @(x) (1+x).^(b+1)./(b+1);
vals_exact = feval(g_exact, x);
err = vals_g - vals_exact;
pass(1) = (norm(err, inf) < 1e1*eps*norm(vals_exact, inf));
    

% fractional pole with order < -1 at the right endpoint:
data.exponents = [0 d];
data.singType = {'none', 'sing'};
f = singfun(@(x) (1-x).^d, data, pref);
g = cumsum(f);
vals_g = feval(g, x);
g_exact = @(x)-(1-x).^(d+1)./(d+1) + 2^(d+1)/(d+1);
vals_exact = feval(g_exact, x);
err = vals_g - vals_exact;
pass(2) = (norm(err, inf) < 1e3*eps*norm(vals_exact, inf));

% fractional root at the left endpoint:
data.exponents = [a 0];
data.singType = {'root', 'none'};
f = singfun(@(x) (1+x).^a, data, pref);
g = cumsum(f);
vals_g = feval(g, x); 
g_exact = @(x) (1+x).^(a+1)./(a+1);
vals_exact = feval(g_exact, x);
err = vals_g - vals_exact;
pass(3) = (norm(err, inf) < eps*norm(vals_exact, inf));

% Integer pole:
data.exponents = [0 -4];
data.singType = {'none', 'pole'};
f = singfun(@(x) (1-x).^(-4), data, pref);
g = cumsum(f);
vals_g = feval(g, x);
g_exact = @(x)(1-x).^(-3)/3 - 2^(-3)/3;
vals_exact = feval(g_exact, x);
err = vals_g - vals_exact;
pass(4) = (norm(err, inf) < 1e2*eps*norm(vals_exact, inf));


%% Tests without closed form solution:

%%
dom = [-1+10^(-D) 1];
scl = (dom(2)-dom(1))/2;

% Fractional pole at the left endpoint:
% Integrate first, then restrict:
data.exponents = [b 0];
data.singType = {'sing', 'none'};
f = singfun(@(x) cos(x.^2+3).*((1+x).^b), data, pref);
u = cumsum(f);
g = restrict(u, dom);

% Restrict first, then integrate:
v = restrict(f, dom);
h = scl*cumsum(v);
h = h - feval(h, 1) + feval(u, 1);

vals_g = feval(g, x);
vals_exact = feval(h, x);
err = norm(vals_g - vals_exact, inf);
tol = 1e4*eps*norm(vals_exact, inf);
pass(5) = (err < tol);

%%
% fractional root at the left endpoint:
dom = [-1+10^(-D) 1];
scl = (dom(2)-dom(1))/2;

% Integrate first, then restrict:
data.exponents = [c 0];
data.singType = {'root', 'none'};
f = singfun(@(x) cos(sin(x)).*(1+x).^c, data, pref);
u = cumsum(f);
g = restrict(u, dom);

% Restrict first, then integrate:
v = restrict(f, dom);
h = scl*cumsum(v);
h = h - feval(h, 1) + feval(u, 1);

vals_g = feval(g, x);
vals_exact = feval(h, x);
err = vals_g - vals_exact;
pass(6) = (norm(err, inf) < 1e2*eps*norm(vals_exact, inf));

% Integer pole at the right endpoint (This can be tested only when log is ready):
% f = singfun(@(x)sin(2+x.^2)./(1-x).^2, [0 -2], {'none', 'pole'}, [], [], pref);
% g = cumsum(f);
% h = cumsum(restrict(f, [-1 1-10^(-D)]));
% vals_g = feval(g, x);
% vals_exact = feval(h, x);
% err = vals_g - vals_exact;
% pass(7) = (norm(err, inf) < eps*norm(vals_exact, inf));

end

% Test file for singfun/cumsum.m

function pass = test_cumsum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
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
e = -2.56;

% Pre-allocate pass matrix
pass = zeros(1, 9);

%%
% Tests with exact solution:

% fractional pole with order > -1 at the left endpoint:
f = singfun(@(x) (1+x).^b, [b 0], {'sing', 'none'}, [], [], pref);
g = cumsum(f);
vals_g = feval(g, x); 
g_exact = @(x) (1+x).^(b+1)./(b+1);
vals_exact = feval(g_exact, x);
err = vals_g - vals_exact;
pass(1) = (norm(err, inf) < get(f,'epslevel')*norm(vals_exact, inf));

% fractional pole with order < -1 at the right endpoint:
f = singfun(@(x) (1-x).^d, [0 d], {'none', 'sing'}, [], [], pref);
g = cumsum(f);
vals_g = feval(g, x);
g_exact = @(x)-(1-x).^(d+1)./(d+1);
vals_exact = feval(g_exact, x);
err = vals_g - vals_exact;
pass(2) = (norm(err, inf) < get(f,'epslevel')*norm(vals_exact, inf));

% fractional root at the left endpoint:
f = singfun(@(x) (1+x).^a, [a 0], {'root', 'none'}, [], [], pref);
g = cumsum(f);
vals_g = feval(g, x); 
g_exact = @(x) (1+x).^(a+1)./(a+1);
vals_exact = feval(g_exact, x);
err = vals_g - vals_exact;
pass(3) = (norm(err, inf) < get(f,'epslevel')*norm(vals_exact, inf));

% fractional pole with order < -1 at the right endpoint with multiple 
% integration:
f = singfun(@(x) (1-x).^(2*d), [0 2*d], {'none', 'sing'}, [], [], pref);
g = cumsum(f, 2);
vals_g = feval(g, x);
g_exact = @(x)(1-x).^(2*(d+1))./((2*d+1)*(2*d+2));
vals_exact = feval(g_exact, x);
err = vals_g - vals_exact;
pass(4) = (norm(err, inf) < get(f,'epslevel')*norm(vals_exact, inf));

% Integer pole:
f = singfun(@(x) (1-x).^(-4), [0 -4], {'none', 'pole'}, [], [], pref);
g = cumsum(f, 1);
vals_g = feval(g, x);
g_exact = @(x)(1-x).^(-3)/3;
vals_exact = feval(g_exact, x);
err = vals_g - vals_exact;
pass(5) = (norm(err, inf) < get(f,'epslevel')*norm(vals_exact, inf));

% Integer root:
f = singfun(@(x) (1+x).^2, [2 0], {'none', 'sing'}, [], [], pref);
g = cumsum(f, 2);
vals_g = feval(g, x);
g_exact = @(x)(1+x).^4/12;
vals_exact = feval(g_exact, x);
err = vals_g - vals_exact;
pass(6) = (norm(err, inf) < get(f,'epslevel')*norm(vals_exact, inf));

%% Tests without closed form solution:

%%
dom = [-1+10^(-D) 1];
scl = (dom(2)-dom(1))/2;

% Fractional pole at the left endpoint:
% Integrate first, then restrict:
f = singfun(@(x) cos(x.^2+3).*((1+x).^b), [b 0], {'sing', 'none'}, [], [], ...
    pref);
u = cumsum(f);
g = restrict(u, dom);

% Restrict first, then integrate:
v = restrict(f, dom);
h = scl*cumsum(v);
h = h - feval(h, 1) + feval(u, 1);

vals_g = feval(g, x);
vals_exact = feval(h, x);
err = vals_g - vals_exact;
pass(7) = (norm(err, inf) < 1e3*get(f,'epslevel')*norm(vals_exact, inf));

%%
% Fractional pole at the right endpoint and multiple integration:
dom = [-1 1-10^(-D)];
scl = (dom(2)-dom(1))/2;

% Integrate first, then restrict:
f = singfun(@(x) sin(x).*((1-x).^e), [0 e], {'none', 'sing'}, [], [], pref);
u = cumsum(f, 2);
C = feval(cumsum(f), -1);  % Get the constant of integration
g = restrict(u, dom);

% Restrict first, then integrate:
v = restrict(f, dom);
vv = scl*cumsum(v) + C;
h = scl*cumsum(vv);
h = h + feval(u, -1);

vals_g = feval(g, x);
vals_exact = feval(h, x);
err = vals_g - vals_exact;
pass(8) = (norm(err, inf) < 1e7*get(f,'epslevel')*norm(vals_exact, inf));

%%
% fractional root at the left endpoint:
dom = [-1+10^(-D) 1];
scl = (dom(2)-dom(1))/2;

% Integrate first, then restrict:
f = singfun(@(x) cos(sin(x)).*(1+x).^c, [c 0], {'root', 'none'}, [], [], pref);
u = cumsum(f);
g = restrict(u, dom);

% Restrict first, then integrate:
v = restrict(f, dom);
h = scl*cumsum(v);
h = h - feval(h, 1) + feval(u, 1);

vals_g = feval(g, x);
vals_exact = feval(h, x);
err = vals_g - vals_exact;
pass(9) = (norm(err, inf) < 1e2*get(f,'epslevel')*norm(vals_exact, inf));

% Integer pole at the right endpoint (This can be tested only when log is ready):
% f = singfun(@(x)sin(2+x.^2)./(1-x).^2, [0 -2], {'none', 'pole'}, [], [], pref);
% g = cumsum(f);
% h = cumsum(restrict(f, [-1 1-10^(-D)]));
% vals_g = feval(g, x);
% vals_exact = feval(h, x);
% err = vals_g - vals_exact;
% pass(10) = (norm(err, inf) < get(f,'epslevel')*norm(vals_exact, inf));

end

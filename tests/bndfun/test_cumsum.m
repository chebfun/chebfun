% Test file for bndfun/cumsum.m

function pass = test_cumsum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

% Set a domain
dom = [-2 7];
a = dom(1);

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

%%
% Spot-check antiderivatives for a couple of functions. We also check that 
% feval(cumsum(f), a) == 0 each time.

f = bndfun(@(x) exp(x/10) - 1, dom, [], [], pref);
F = cumsum(f);
F_ex = @(x) 10*exp(x/10) - x;
err = feval(F, x) - F_ex(x);
pass(1) = (norm(diff(err), inf) < 100*get(f, 'vscale')*get(f, 'epslevel')) ...
    && (abs(feval(F, a)) <= get(f, 'vscale')*get(f, 'epslevel'));

f = bndfun(@(x) 1./(1 + x.^2), dom, [], [], pref);
F = cumsum(f);
F_ex = @(x) atan(x);
err = feval(F, x) - F_ex(x);
pass(2) = (norm(diff(err), inf) < 100*get(f, 'vscale')*get(f, 'epslevel')) ...
    && (abs(feval(F, a)) <= get(f, 'vscale')*get(f, 'epslevel'));

f = bndfun(@(x) cos(1e4*x), dom, [], [], pref);
F = cumsum(f);
F_ex = @(x) sin(1e4*x)/1e4;
err = feval(F, x) - F_ex(x);
pass(3) = (norm(diff(err), inf) < 100*get(f, 'vscale')*get(f, 'epslevel')) ...
    && (abs(feval(F, a)) <= get(f, 'vscale')*get(f, 'epslevel'));

z = exp(2*pi*1i/6);
f = bndfun(@(t) sinh(t*z), dom, [], [], pref);
F = cumsum(f);
F_ex = @(t) cosh(t*z)/z;
err = feval(F, x) - F_ex(x);
pass(4) = (norm(diff(err), inf) < 100*get(f, 'vscale')*get(f, 'epslevel')) ...
    && (abs(feval(F, a)) <= get(f, 'vscale')*get(f, 'epslevel'));

%%
% Check that applying cumsum() and direct construction of the antiderivative
% give the same results (up to a constant).

f = bndfun(@(x) sin(4*x).^2, dom, [], [], pref);
F = bndfun(@(x) 0.5*x - 0.0625*sin(8*x), dom, [], [], pref);
G = cumsum(f);
err = feval(G - F, x);
pass(5) = (norm(diff(err), inf) < 10*get(f, 'vscale')*get(f, 'epslevel')) && ...
    (abs(feval(G, a)) < 10*get(f, 'vscale')*get(f, 'epslevel'));

%%
% Check that diff(cumsum(f)) == f and that cumsum(diff(f)) == f up to a
% constant.

f = bndfun(@(x) x.*(x - 1).*sin(x), dom, [], [], pref);
g = diff(cumsum(f));
tol_f = 10*get(f, 'vscale')*get(f, 'epslevel');
tol_g = 10*get(g, 'vscale')*get(g, 'epslevel');

err = feval(f, x) - feval(g, x);
pass(6) = (norm(diff(err), inf) < 10*max(tol_f, tol_g));
h = cumsum(diff(f));
err = feval(f, x) - feval(h, x);
pass(7) = (norm(diff(err), inf) < 10*max(tol_f, tol_g) && ...
    (abs(feval(h, a)) < 10*max(tol_f, tol_g)));

%%
% Check operation for array-valued bndfun objects.

f = bndfun(@(x) [sin(x) x.^2 exp(1i*x)], dom, [], [], pref);
F_exact = bndfun(@(x) [-cos(x) x.^3/3 exp(1i*x)/1i], dom, [], [], pref);
F = cumsum(f);
err = feval(F, x) - feval(F_exact, x);
pass(8) = (norm(diff(err), inf) < ...
    10*max(get(f, 'vscale').*get(f, 'epslevel'))) && ...
    all(abs(feval(F, a)) < max(get(f, 'vscale').*get(f, 'epslevel')));

%% Test on singular function:

% Singularity at one endpoint:

dom = [-2 7];
pow = -0.64;
op = @(x) (x-dom(1)).^pow;
pref.singPrefs.exponents = [pow 0];
f = bndfun(op, dom, [], [], pref);
g = cumsum(f);
vals_g = feval(g, x); 
g_exact = @(x) (x-a).^(pow+1)./(pow+1);
vals_exact = feval(g_exact, x);
err = vals_g - vals_exact;
pass(9) = ( norm(err, inf) < 10*get(f,'epslevel')*norm(vals_exact, inf) );

% Singularities at both endpoints:
mid = mean(dom);
x = diff(dom) * rand(100, 1) + dom(1);
pow1 = -0.5;
pow2 = -0.5;
op = @(x) (x-dom(1)).^pow1.*(dom(2)-x).^pow2;
pref.singPrefs.exponents = [pow1 pow2];
f = bndfun(op, dom, [], [], pref);

% We temporarily disable this warning: 
warning('off', 'CHEBFUN:SINGFUN:plus');
g = cumsum(f);
warning('on', 'CHEBFUN:SINGFUN:plus');
x = sort(x);
x1 = x(x <= mid);
x2 = x(x > mid);
vals_g1 = feval(g{1}, x1); 
vals_g2 = feval(g{2}, x2);

% This exact solution is acquired using Matlab symbolic toolbox:
g_exact = @(x) 4*atan((sqrt(dom(2)) - sqrt(dom(2) - x))./ ...
    (sqrt(x - dom(1)) - sqrt(-dom(1)))); 
g_exact1 = @(x) g_exact(x) - g_exact(dom(1)) + get(g{1}, 'lval');
g_exact2 = @(x) g_exact(x) - g_exact(dom(2)) + get(g{2}, 'rval');
vals_exact1 = feval(g_exact1, x1);
vals_exact2 = feval(g_exact2, x2);
err1 = vals_g1 - vals_exact1;
err2 = vals_g2 - vals_exact2;
pass(10) = ( norm(err1, inf) < 1e3*get(g{1},'epslevel')* ...
    norm(vals_exact1, inf) ) && ...
    ( norm(err2, inf) < 1e3*get(g{2},'epslevel')*norm(vals_exact2, inf) );

end
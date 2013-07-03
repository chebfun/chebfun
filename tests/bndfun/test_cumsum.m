% Test file for bndfun/cumsum.m

function pass = test_cumsum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = fun.pref;
end

% Set a tolerance.
tol = pref.fun.eps;

% Set a domain
dom = [-2 7];
a = dom(1);

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

pass = zeros(1, 10); % Pre-allocate pass matrix

%%
% Spot-check antiderivatives for a couple of functions. We also check that 
% feval(cumsum(f), a) == 0 each time.

f = bndfun(@(x) exp(x/10) - 1, dom, [], [], pref);
F = cumsum(f);
F_ex = @(x) 10*exp(x/10) - x;
err = feval(F, x) - F_ex(x);
pass(1) = (norm(diff(err), inf) < 25*get(f, 'vscale')*tol) && ...
    (abs(feval(F, a)) < get(f, 'vscale')*tol);

f = bndfun(@(x) 1./(1 + x.^2), dom, [], [], pref);
F = cumsum(f);
F_ex = @(x) atan(x);
err = feval(F, x) - F_ex(x);
pass(2) = (norm(diff(err), inf) < 30*get(f, 'vscale')*tol) && ...
    (abs(feval(F, a)) < get(f, 'vscale')*tol);

f = bndfun(@(x) cos(1e4*x), dom, [], [], pref);
F = cumsum(f);
F_ex = @(x) sin(1e4*x)/1e4;
err = feval(F, x) - F_ex(x);
pass(3) = (norm(diff(err), inf) < 1e4*get(f, 'vscale')*tol) && ...
    (abs(feval(F, a)) < get(f, 'vscale')*tol);

z = exp(2*pi*1i/6);
f = bndfun(@(t) sinh(t*z), dom, [], [], pref);
F = cumsum(f);
F_ex = @(t) cosh(t*z)/z;
err = feval(F, x) - F_ex(x);
pass(4) = (norm(diff(err), inf) < 20*get(f, 'vscale')*tol) && ...
    (abs(feval(F, a)) < get(f, 'vscale')*tol);

%%
% Check that applying cumsum() and direct construction of the antiderivative
% give the same results (up to a constant).

f = bndfun(@(x) sin(4*x).^2, dom, [], [], pref);
F = bndfun(@(x) 0.5*x - 0.0625*sin(8*x), dom, [], [], pref);
G = cumsum(f);
err = feval(G - F, x);
pass(5) = (norm(diff(err), inf) < 20*get(f, 'vscale')*tol) && ...
    (abs(feval(G, a)) < get(f, 'vscale')*tol);

%%
% Check that diff(cumsum(f)) == f and that cumsum(diff(f)) == f up to a
% constant.

f = bndfun(@(x) x.*(x - 1).*sin(x), dom, [], [], pref);
g = diff(cumsum(f));
err = feval(f, x) - feval(g, x);
pass(6) = (norm(diff(err), inf) < 20*get(g, 'vscale')*get(f, 'vscale')*tol);
h = cumsum(diff(f));
err = feval(f, x) - feval(h, x);
pass(7) = (norm(diff(err), inf) < 30*get(g, 'vscale')*get(f, 'vscale')*tol) && ...
    (abs(feval(h, a)) < get(g, 'vscale')*get(f, 'vscale')*tol);

%%
% Check operation for array-valued bndfun objects.

f = bndfun(@(x) [sin(x) x.^2 exp(1i*x)], dom, [], [], pref);
F_exact = bndfun(@(x) [-cos(x) x.^3/3 exp(1i*x)/1i], dom, [], [], pref);
F = cumsum(f);
err = feval(F, x) - feval(F_exact, x);
pass(8) = (norm(diff(err), inf) < max(get(f, 'vscale'))*get(f, 'epslevel')) && ...
    all(abs(feval(F, a)) < max(get(f, 'vscale'))*get(f, 'epslevel'));

%%
% Check operation for second and third order cumsums.
f = bndfun(@(x) sin(x), dom, [], [], pref);
F2_exact = bndfun(@(x) -sin(x) + x.*cos(a), dom, [], [], pref);

F2 = cumsum(f, 2);
err = feval(F2, x) - feval(F2_exact, x);
pass(9) = (norm(diff(err), inf) < 2*get(F2, 'vscale')^2*get(f, 'epslevel')) && ...
    abs(feval(F2, a)) < get(F2, 'vscale')^2*get(f, 'epslevel');

F3_exact = bndfun(@(x) cos(x) + x.^2*cos(a)/2 + x*(sin(a) - a*cos(a)), ...
    dom, [], [], pref);
F3 = cumsum(f, 3);
err = feval(F3, x) - feval(F3_exact, x);
pass(10) = (norm(diff(err), inf) < get(F3, 'vscale')^3*get(f, 'epslevel')) && ...
    abs(feval(F3, a)) < get(F3, 'vscale')^3*get(f, 'epslevel');

end

% Test file for bndfun/cumsum.m

function pass = test_cumsum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = bndfun.pref;
end
pref = chebtech.pref(pref);

% Set a tolerance.
tol = pref.bndfun.eps;

% Set a domain
dom = [-2 7];
a = dom(1);
b = dom(2);

% Generate a few random points to use as test values.
seedRNG(6178);
x = diff(dom) * rand(100, 1) + dom(1);

pass = zeros(1, 10); % Pre-allocate pass matrix

%%
% Spot-check antiderivatives for a couple of functions. We also check that 
% feval(cumsum(f), a) == 0 each time.

f = bndfun(@(x) exp(x/10) - 1, dom, [], [], pref);
F = cumsum(f);
F_ex = @(x) 10*exp(x/10) - x - (10*exp(a/10) - a);
err = feval(F, x) - F_ex(x);
pass(1) = (norm(err, inf) < 5*f.onefun.vscale*tol) && ...
    (abs(feval(F, a)) < f.onefun.vscale*tol);

f = bndfun(@(x) 1./(1 + x.^2), dom, [], [], pref);
F = cumsum(f);
F_ex = @(x) atan(x)-atan(a);
err = feval(F, x) - F_ex(x);
pass(2) = (norm(err, inf) < 20*f.onefun.vscale*tol) && ...
    (abs(feval(F, a)) < f.onefun.vscale*tol);

f = bndfun(@(x) cos(1e4*x), dom, [], [], pref);
F = cumsum(f);
F_ex = @(x) (sin(1e4*x)-sin(1e4*a))/1e4;
err = feval(F, x) - F_ex(x);
pass(3) = (norm(err, inf) < 1e4*f.onefun.vscale*tol) && ...
    (abs(feval(F, a)) < f.onefun.vscale*tol);

z = exp(2*pi*1i/6);
f = bndfun(@(t) sinh(t*z), dom, [], [], pref);
F = cumsum(f);
F_ex = @(t) cosh(t*z)/z - (cosh(a*z)/z);
err = feval(F, x) - F_ex(x);
pass(4) = (norm(err, inf) < 20*f.onefun.vscale*tol) && ...
    (abs(feval(F, a)) < f.onefun.vscale*tol);

%%
% Check that applying cumsum() and direct construction of the antiderivative
% give the same results (up to a constant).

f = bndfun(@(x) sin(4*x).^2, dom, [], [], pref);
F = bndfun(@(x) 0.5*x - 0.0625*sin(8*x) - (a/2-sin(8*a)/16), dom, [], [], pref);
G = cumsum(f);
err = feval(G - F, x);
pass(5) = (norm(err, inf) < 3*f.onefun.vscale*tol) && ...
    (abs(feval(G, a)) < f.onefun.vscale*tol);

%%
% Check that diff(cumsum(f)) == f and that cumsum(diff(f)) == f up to a
% constant.

f = bndfun(@(x) x.*(x - 1).*sin(x) - (a*(a-1)*sin(a)), dom, [], [], pref);
g = diff(cumsum(f));
err = feval(f, x) - feval(g, x);
pass(6) = (norm(err, inf) < 20*g.onefun.vscale*f.onefun.vscale*tol);
h = cumsum(diff(f));
err = feval(f, x) - feval(h, x);
pass(7) = (norm(err, inf) < 30*g.onefun.vscale*f.onefun.vscale*tol) && ...
    (abs(feval(h, a)) < g.onefun.vscale*f.onefun.vscale*tol);

%%
% Check operation for array-valued bndfun objects.

f = bndfun(@(x) [sin(x) x.^2 exp(1i*x)], dom, [], [], pref);
F_exact = bndfun(@(x) [(-cos(x)+cos(a)) (x.^3/3-a^3/3) ...
    ((exp(1i*x)-exp(1i*a))/1i)], dom, [], [], pref);
F = cumsum(f);
err = feval(F, x) - feval(F_exact, x);
pass(8) = (norm(err, inf) < 5*max(f.onefun.vscale)*tol) && ...
    all(abs(feval(F, a) < 2*max(f.onefun.vscale)*tol));

%%
% Check operation for second and third order cumsums.
f = bndfun(@(x) sin(x), dom, [], [], pref);
F2_exact = bndfun(@(x) -sin(x)+x.*cos(a)+sin(a)-a*cos(a), ...
    dom, [], [], pref);
F2 = cumsum(f,2);
err = feval(F2, x) - feval(F2_exact, x);
pass(9) = (norm(err, inf) < 5*(f.onefun.vscale)^2*tol) && ...
    abs(feval(F2, a) < (f.onefun.vscale)^2*tol);

F3_exact = bndfun(@(x) cos(x)+x.^2*cos(a)/2+x*(sin(a)-a*cos(a)) + ...
    (-cos(a)+a^2*cos(a)/2-a*sin(a)), dom, [], [], pref);
F3 = cumsum(f,3);
err = feval(F3, x) - feval(F3_exact, x);
pass(10) = (norm(err, inf) < 2*(f.onefun.vscale)^3*tol) && ...
    abs(feval(F3, a) < (f.onefun.vscale)^3*tol);

end
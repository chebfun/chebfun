% Test file for bndfun/cumsum.m

function pass = test_cumsum(pref)

% Get preferences.
if ( nargin < 1 )
    pref = bndfun.pref;
end
pref = chebtech.pref(pref);

% Set a tolerance.
tol = 50*pref.chebtech.eps;

% Set a domain
dom = [-2 7];

% Generate a few random points to use as test values.
seedRNG(6178);
x = (diff(dom)/2) * rand(100, 1) + mean(dom);

pass = zeros(1, 10); % Pre-allocate pass matrix

%%
% Spot-check antiderivatives for a couple of functions.  We verify that the
% bndfun antiderivatives match the true ones up to a constant by checking
% that the standard deviation of the difference between the two on a large
% random grid is small. We also check that feval(cumsum(f), -1) == 0 each
% time.

f = bndfun(@(x) exp(x/10) - 1, dom, [], [], pref);
F = cumsum(f);
F_ex = @(x) 10*exp(x/10) - x;
err = feval(F, x) - F_ex(x);
pass(n, 1) = (std(err) < tol) && (abs(feval(F, -1)) < tol);

f = bndfun(@(x) 1./(1 + x.^2), dom, [], [], pref);
F = cumsum(f);
F_ex = @(x) atan(x);
err = feval(F, x) - F_ex(x);
pass(n, 2) = (std(err) < tol) && (abs(feval(F, -1)) < tol);

f = bndfun(@(x) cos(1e4*x), dom, [], [], pref);
F = cumsum(f);
F_ex = @(x) sin(1e4*x)/1e4;
err = feval(F, x) - F_ex(x);
pass(n, 3) = (std(err) < tol) && (abs(feval(F, -1)) < tol);

z = exp(2*pi*1i/6);
f = bndfun(@(t) sinh(t*z), dom, [], [], pref);
F = cumsum(f);
F_ex = @(t) cosh(t*z)/z;
err = feval(F, x) - F_ex(x);
pass(n, 4) = (std(err) < tol) && (abs(feval(F, -1)) < tol);

%%
% Check that applying cumsum() and direct construction of the antiderivative
% give the same results (up to a constant).

f = bndfun(@(x) sin(4*x).^2, dom, [], [], pref);
F = bndfun(@(x) 0.5*x - 0.0625*sin(8*x), dom, [], [], pref);
G = cumsum(f);
err = G - F;
pass(n, 5) = (std(err.onefun.values) < tol) && (abs(feval(G, -1)) < tol);

%%
% Check that diff(cumsum(f)) == f and that cumsum(diff(f)) == f up to a
% constant.

f = bndfun(@(x) x.*(x - 1).*sin(x) + 1, dom, [], [], pref);
g = diff(cumsum(f));
err = feval(f, x) - feval(g, x);
pass(n, 6) = (norm(err, inf) < 100*tol);
h = cumsum(diff(f));
err = feval(f, x) - feval(h, x);
pass(n, 7) = (std(err) < tol)  && (abs(feval(h, -1)) < tol);

%%
% Check operation for array-valued chebtech objects.

f = bndfun(@(x) [sin(x) x.^2 exp(1i*x)], dom, [], [], pref);
F_exact = bndfun(@(x) [(-cos(x)) (x.^3/3) (exp(1i*x)/1i)], dom, [], [], pref);
F = cumsum(f);
err = std(feval(F, x) - feval(F_exact, x));
pass(n, 8) = (norm(err, inf) < tol)  && all(abs(feval(F, -1)) < tol);

%%
% Check operation for second and third order cumsums.
f = bndfun(@(x) sin(x), dom, [], [], pref);
F2_exact = bndfun(@(x) -sin(x)+x.*cos(dom(1))+sin(dom(1))+cos(dom(1)), ...
    dom, [], [], pref);
F2 = cumsum(f,2);
err = std(feval(F2, x) - feval(F2_exact, x));
pass(n, 9) = (norm(err, inf) < tol)  && abs(feval(F2, dom(1)) < tol);

F3_exact = bndfun(@(x) cos(x)+x.^2*cos(dom(1))/2+x*(sin(dom(1))+cos(dom(1))) - ...
    (cos(dom(1))/2-sin(dom(1))), [], [], pref);
F3 = cumsum(f,3);
err = std(feval(F3, x) - feval(F3_exact, x));
pass(n, 10) = (norm(err, inf) < tol)  && abs(feval(F3, dom(1)) < tol);

end
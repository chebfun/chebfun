function pass = test_wronskian(pref)

% Choose a tolerance:
tol = 1e-11;

% Check different input styles for a simple problem:
op = @(x,u) diff(u, 2) + (u);
L = chebop(op, [-pi, pi]);
x = chebfun('x', [-pi, pi]);
w = wronskian(L, cos(x), sin(x));
pass(1) = norm(w - 1, inf) < tol;
w = wronskian(L, [cos(x), sin(x)]);
pass(2) = norm(w - 1, inf) < tol;
w = wronskian(L, chebmatrix({cos(x), sin(x)}));
pass(3) = norm(w - 1, inf) < tol;


% A slightly harder problem:
op = @(u) diff(u,2) + 4*diff(u) - 5*u;
dom = [0, 5];
L = chebop(op, dom);
x = chebfun(@(x) x, dom);
g = exp(x);
f = exp(-5*x);
w = wronskian(L, f, g);
% Compare with the exact wronskian:
pass(4) = norm(w - 6*exp(-4*x), inf) < tol;

end

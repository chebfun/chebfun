function pass = test_zerothOrder(~)
% TEST_ZEROTHORDER   Solve problems without any differentiation/integration

%% Do a manual Newton iteration to check everything is working at lower levels
x = chebfun('x');
% RHS
f = sin(x)+2;
% Initial guess
u = 0*x;
op = @(u) u.^2 + sin(u) + exp(u) - f;
p = cheboppref();
p.discretization = @chebcolloc2;
normdel = inf;
while normdel > 1e-10
    uad = adchebfun(u);
    res = op(uad);
    del = linsolve(linop(res.jacobian), res.func, p);
    normdel = norm(del);
    u = chebfun(u);
    u = u - del;
end

pass(1) = normdel < 1e-10;

%% Iteration with CHEBOPs
x = chebfun('x');
f = sin(4*x) + 2;
N = chebop(@(u) 4*u.^2 - x+.1+sin(u));
u = N\f;
pass(2) = norm(N(u) - f) < 1e-10;

%% Coupled systems
x = chebfun('x');
f = sin(4*x) + 2;
N = chebop(@(x,u,v) [u.^2 - x+.1+sin(v); v + exp(v) + u]);
[u, v] = N\[f; f];
pass(3) = norm(N(u,v) - [f;f]) < 1e-10;

end
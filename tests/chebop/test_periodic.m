function pass = test_periodic(pref)
% Test 'periodic' syntax for CHEBOP

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-10;

%%

% A simple ODE:
N = chebop(@(u) diff(u, 2) - u);
% Periodic RHS:
rhs = chebfun(@(x) sin(pi*x));

% Apply periodic BCs manually:
N.bc = @(u) [u(-1) - u(1) ; feval(diff(u),-1) - feval(diff(u), 1)];
u = mldivide(N, rhs, pref);

% Apply via 'periodic':
N.bc = 'periodic';
v = mldivide(N, rhs, pref);

% Compare results:
err(1) = norm(u - v);

%%

% A periodic piecewise system:
d = [-pi, 0, pi];
A = chebop(d);
A.op = @(x, u, v) [u-diff(v) ; diff(u,2)+v];
x = chebfun('x',d);
f = [chebfun(0, d) ; cos(x)];
A.bc = 'periodic';
uv = mldivide(A, f, pref);
u = uv{1}; v = uv{2};

trueSoln = [cos(x+3*pi/4), cos(x+pi/4)]/sqrt(2);
err(2) = norm([u v] - trueSoln);

%%

pass = err < tol;

end


function pass = test_periodic(pref)
% Test 'periodic' syntax for CHEBOP

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-10;

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
err = norm(u - v);
pass = err < tol;

end


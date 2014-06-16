function pass = test_bc(pref)
% Solve two simple linear BVPs, check the residual and the precision of the
% boundary conditions.

% Note: These tests are taken from chebop_bc and chebop_bc2 in V4.

if ( nargin == 0 )
    pref = cheboppref;
end
pref.errTol = 1e-12;

%% linear: numeric and string input
d = [-3 4];
A = chebop(@(x, u) diff(u,2) + 4*diff(u) + u, d);
A.lbc = -1;
A.rbc = 'neumann';
f = chebfun( 'exp(sin(x))', d );
u = solvebvp(A, f, pref);

err(1) = norm(diff(u,2) + 4*diff(u) + u - f);
err(2) = abs(u(d(1)) + 1);
err(3) = abs(feval(diff(u), d(2)));

%% Linear: function_handle input
d = [-1 0];
A = chebop(@(x,u) diff(u,2) + 4*diff(u) + 200*u, d);
A.lbc = @(u) [diff(u)+2*u-1];
A.rbc = @(u) diff(u);
f = chebfun( 'x.*sin(3*x).^2',d );
u = solvebvp(A, f, pref);
du = diff(u);

err(4) = norm(diff(u,2) + 4*diff(u) + 200*u - f);
err(5) = abs(du(d(1)) + 2*u(d(1)) - 1);
err(6) = abs(feval(diff(u), d(2)));

%% Linear: exotic bcs:
dom = [-2 1];
N = chebop(@(x,u) diff(u,2) + u, dom);
N.bc = @(x,u) [feval(diff(u),0) ; sum(u)];
x = chebfun(@(x) x, dom);
rhs = sin(x);
u = solvebvp(N, rhs, pref);

err(7) = norm( N(x,u) - rhs );
err(8) = norm(N.bc(x,u));

%% Nonlinear

dom = [-1 1];
N = chebop(@(x,u) diff(u,2) + sin(u), dom);
N.bc = @(x,u) [feval(diff(u),0) ; sum(u)];
x = chebfun(@(x) x, dom);
rhs = sin(x);
u = solvebvp(N, rhs, pref);

err(9) = norm( N(x,u) - rhs );
err(10) = norm(N.bc(x, u));

%%

tol = 1e-10;
pass = err < tol;

end

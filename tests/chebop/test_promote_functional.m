function pass = test_promote_functional(~)

% Test operator + functional and operator*functional.

tol = 1e-10;

%%
N = chebop(@(u) diff(u,2) + sum(u));
N.lbc = 0;
N.rbc = 0;
u = N\1;
pass(1) = norm(N(u) - 1, inf) < tol;

%%
N = chebop(@(u) diff(u,2) + u.*sum(u));
N.lbc = 0;
N.rbc = 0;
u = N\1;
pass(2) = norm(N(u) - 1, inf) < tol;

%%
N = chebop(@(x,u,v) [diff(u) ; diff(v)]);
N.bc = @(x,u, v) [u(0) - 1; u(0).*v(.5)];
N.init = [1 ; 0];
pass(3) = size(matrix(linearize(N), 3), 1) == 8;
u = N\1;
pass(4) = norm(N(u) - 1, inf) < tol;

%%

N = chebop(@(x,u,v) [diff(u) ; diff(v) + v.*sum(u) + u(.3)]);
N.bc = @(x,u, v) [u(0)-1; u(0).*v(.5)];
pass(5) = size(matrix(linearize(N), 3), 1) == 8;
u = N\1;
pass(6) = norm(N(u) - 1, inf) < tol;

%% Periodic (with values):

N = chebop(@(u) diff(u,2) + sum(u));
N.bc = 'periodic';
rhs = chebfun(@(x) cos(pi*x));
u = N\rhs;
pass(7) = norm(N(u) - rhs, inf) < tol;
pass(8) = isequal(get(u.funs{1}, 'tech'), @trigtech);

%% Periodic (with coeffs):

N = chebop(@(u) diff(u,2) + sum(u));
N.bc = 'periodic';
rhs = chebfun(@(x) cos(pi*x));
pref = cheboppref();
pref.discretization = 'coeffs';
u = solvebvp(N, rhs, pref);
pass(9) = norm(N(u) - rhs, inf) < tol;
pass(10) = isequal(get(u.funs{1}, 'tech'), @trigtech);

end

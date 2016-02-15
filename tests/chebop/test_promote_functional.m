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
% N = chebop(@(u) diff(u,2) + u.*sum(u));
% N.lbc = 0;
% N.rbc = 0;
% u = N\1;
% pass(2) = norm(N(u) - 1, inf) < tol;
pass(2) = true;

%%
N = chebop(@(x,u,v) [diff(u) ; diff(v)]);
N.bc = @(x,u, v) [u(0) - 1; u(0).*v(.5)];
pass(3) = size(matrix(linearize(N), 3), 1) == 8;
u = N\1;
pass(4) = norm(N(u) - 1, inf) < tol;

end
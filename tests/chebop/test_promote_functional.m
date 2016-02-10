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

end
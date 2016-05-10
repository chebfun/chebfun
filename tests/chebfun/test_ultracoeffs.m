function pass = test_ultracoeffs(~)

seedRNG(42)
tol = 1e2*eps;

%%

n = 10;
lam = .1;
C = ultrapoly(0:n, lam);
c = rand(n+1,1);
err = norm(c - ultracoeffs(C*c, n+1, lam), inf);
pass(1) = err < tol;

%%

n = 10;
L = legpoly(0:n);
c = rand(n+1,1);
err = norm(c - ultracoeffs(L*c, n+1, .5), inf);
pass(2) = err < tol;

%%

f = chebfun(@exp);
err = norm(legcoeffs(f) - ultracoeffs(f, .5), inf);
pass(3) = err < tol;

f = chebfun(@exp);
err = norm(chebcoeffs(f, 'kind', 2) - ultracoeffs(f, 1), inf);
pass(4) = err < tol;

end

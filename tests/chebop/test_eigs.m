function pass = test_eigs(pref)

% NOTE: This was taken chebop_eigs in the V4 tests.

if ( nargin == 0 )
    pref = cheboppref();
end

%% With linops
d = domain(0, pi);
L = -diff(d, 2);
e = functionalBlock.eval(d);
L = addConstraint(L,chebmatrix(e(0)));
L = addConstraint(L,chebmatrix(e(pi)));
[~, D] = eigs(L, 10, pref);
e1 = sqrt(-diag(D));
err(1) = norm(e1 - 1i*(1:10).',inf);

%% With chebops
d = [0, pi];
N = chebop(d);
N.op = @(x,u) -diff(u,2);
N.lbc = 'dirichlet';
N.rbc = 'dirichlet';
[~, D] = eigs(N, 10, pref);
e2 = sqrt(-diag(D));
err(2) = norm(e2 - 1i*(1:10).',inf);

pass = err < 1e-10;

end
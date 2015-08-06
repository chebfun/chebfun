function pass = test_eigs_basic(pref)

% NOTE: This was taken chebop_eigs in the V4 tests.

if ( nargin == 0 )
    pref = cheboppref();
end

tolVals = 1e-10;
tolFuns = 1e-7;
%% With linops
d = domain(0, pi);
L = -diff(d, 2);
e = functionalBlock.eval(d);
L = addConstraint(L,chebmatrix(e(0)));
L = addConstraint(L,chebmatrix(e(pi)));
pref.discretization = @chebcolloc2;
[V, D] = eigs(L, 10, pref);
e1 = sqrt(-diag(D));
errVals(1) = norm(e1 - 1i*(1:10).',inf);
errFuns(1) = norm(L*V-V*D);
%% With chebops
d = [0, pi];
N = chebop(d);
N.op = @(x,u) -diff(u,2);
N.lbc = 'dirichlet';
N.rbc = 'dirichlet';
pref.discretization = 'values';
[V, D] = eigs(N, 10, pref);
e2 = sqrt(-diag(D));
errVals(2) = norm(e2 - 1i*(1:10).',inf);
errFuns(2) = norm(N(V)-V*D);
%%

pass = [(errVals < tolVals), (errFuns < tolFuns)];

end

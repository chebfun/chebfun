function pass = test_sparse

% Test sparse collocation discretization on the same example as test_system3.

% Store preference state:
savedPrefs = cheboppref;

% Force sparse discretization:
cheboppref.setDefaults('sparse', true);
pref = cheboppref;

N = chebop(@(x,u,v) [diff(u,2) + sin(v) ; cos(u) + diff(v,2)], [-1:1/2:1]);
N.lbc = @(u,v) [u - 2 ; v - 1]; 
N.rbc = @(u,v) [u - 2 ; v + 1];

tic
[u, v, info1] = solvebvp(N, [0 ; 0], pref);
pass(1) = abs(feval(v, .2) - -0.371250985730553) < 1e-8;
time_sparse = toc;

J = matrix(linearize(N, [u;v]),10);
pass(2) = issparse(J);

% Compare to nonsparse discretisation
cheboppref.setDefaults('sparse', false);
pref = cheboppref;

tic
[u2, v2, info1] = solvebvp(N, [0 ; 0], pref);
time_dense = toc;

pass(3) = norm(u2-u,inf) + norm(v2-v,inf) < 1e-10;

% Revert preference state:
cheboppref.setDefaults(savedPrefs);

% [time_sparse, time_dense]

end
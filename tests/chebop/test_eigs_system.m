function pass = test_eigs_system(pref)
% Eigenvalue test, inspired by Maxwell's equation. The eigenvalues are
% computed first on a global domain, then on a piecewise domain.

% NOTE: This was taken chebop_systemeig in the V4 tests.

if ( nargin == 0 )
    pref = cheboppref();
end

lamCorrect = [ 0
               -.5 + sqrt(3)/2*1i
               -.5 - sqrt(3)/2*1i
               -.5 + sqrt(15)/2*1i
               -.5 - sqrt(15)/2*1i];

%% Smooth domain
d = [0, pi];
A = chebop(@(x, u, v) [-u + diff(v) ; diff(u)], d);
A.lbc = @(u, v) u;
A.rbc = @(u, v) u;
[~, D] = eigs(A, 5, pref);
lam = diag(D);
err(1) = norm( lam - lamCorrect, inf );

%% Piecewise domain
A.domain = [0, pi/2, pi];
[~, D] = eigs(A, 5, pref);
lam_pw = diag(D);
err(2) = norm( lam_pw - lamCorrect, inf )

%%

pass = err < 1e-12;

end
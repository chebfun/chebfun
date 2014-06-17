function pass = test_eigs_system(pref)
% Eigenvalue test, inspired by Maxwell's equation. The eigenvalues are
% computed first on a global domain, then on a piecewise domain.

% NOTE: This was taken chebop_systemeig in the V4 tests.

% NOTE: We compare ABS values to avoid confusions over ordering of complex
% values.

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
[V, D] = eigs(A, 5, pref);
lam = diag(D);
lam = abs(lam);
lamCorrect = abs(lamCorrect);
AV = [];
for vCounter = 1:size(V, 2)
    AV = [AV, A(V(:, vCounter))];
end
err(1) = norm( lam - lamCorrect, inf );
err(2) = norm(AV - V*D);
%% Piecewise domain
A.domain = [0, pi/2, pi];
[V, D] = eigs(A, 5, pref);
lam_pw = diag(D);
lam_pw = abs(lam_pw);
AV = [];
for vCounter = 1:size(V, 2)
    AV = [AV, A(V(:, vCounter))];
end
err(3) = norm( lam_pw - lamCorrect, inf );
err(4) = norm(AV - V*D);
%%

pass = err < 1e-11;

end

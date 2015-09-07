% Test file for chebtech/chebTcoeffs2chebUcoeffs.m.

function pass = test_chebTcoeffs2chebUcoeffs(pref)

% Preference input does not matter.

% Check empty input.
pass(1) = isempty(chebtech.chebTcoeffs2chebUcoeffs([]));

% Check a column vector input.
cT = [1.875 ; 1.75 ; 1.0 ; 0.25 ; 0.125];   % 1 + x + x^2 + x^3 + x^4.
cU = chebtech.chebTcoeffs2chebUcoeffs(cT);
cU_exact = [1.375 ; 0.75 ; 0.4375 ; 0.125 ; 0.0625];
err = norm(cU - cU_exact, Inf);
tol = 10*eps;
pass(2) = err < tol;

% Check a matrix input.
cT = eye(5, 5);
cU = chebtech.chebTcoeffs2chebUcoeffs(cT);  % Generates recurrence matrix.
cU_exact = diag(0.5*ones(5, 1)) + diag(-0.5*ones(3, 1), 2);
cU_exact(1,1) = 1;
err = norm(cU(:) - cU_exact(:), Inf);
tol = 10*eps;
pass(3) = err < tol;

end

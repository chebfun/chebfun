function pass = test_eigs()
% TAD, 10 Jan 2014

tol = 1e-8;

dom = [-pi/2, pi/2];
D2 = operatorBlock.diff(dom, 2);
E = functionalBlock.eval(dom);
El = E(dom(1));
Er = E(dom(end));
L = linop(D2);
L = addbc(L, El, 0);
L = addbc(L, Er, 0);

e_true = flipud(-(1:6).'.^2);

%%
prefs = cheboppref;
prefs.discretization = @chebcolloc2;
[V, D] = eigs(L, 6, prefs);
e = diag(D);
err(1) = norm(e - e_true, inf);
% Check that we actually computed eigenfunctions
err(2) = norm(L*V-V*D);
%%
prefs.discretization = @ultraS;
[V, D] = eigs(L, 6, 0, prefs);
e = diag(D);
err(3) = norm(e - e_true, inf);
err(4) = norm(L*V-V*D);
%%
prefs.discretization = @chebcolloc1;
[V, D] = eigs(L, 6, prefs);
e = diag(D);
err(5) = norm(e - e_true, inf);
err(6) = norm(L*V-V*D);
%%
pass = err < tol;

end

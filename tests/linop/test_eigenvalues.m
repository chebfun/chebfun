function pass = test_eigenvalues()
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
prefs.discretization = @colloc2;
[V, D] = eigs(L, 6, prefs);
e = diag(D);
err(1) = norm(e - e_true, inf);
% Check that we actually computed eigenfunctions
LV = L*V;
err(2) = norm(LV{1}-V{1}*D);
%%
prefs.discretization = @ultraS;
[V, D] = eigs(L, 6, 0, prefs);
e = diag(D);
err(3) = norm(e - e_true, inf);
LV = L*V;
err(4) = norm(LV{1}-V{1}*D);
%%
prefs.discretization = @colloc1;
[V, D] = eigs(L, 6, prefs);
e = diag(D);
err(5) = norm(e - e_true, inf);
LV = L*V;
err(6) = norm(LV{1}-V{1}*D);
%%
pass = err < tol;

end
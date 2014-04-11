function pass = test_eigenvalues()
% TAD, 10 Jan 2014

dom = [-pi/2, pi/2];
D2 = operatorBlock.diff(dom, 2);
E = functionalBlock.eval(dom);
El = E(dom(1));
Er = E(dom(end));
L = linop(D2);
L = addbc(L, El, 0);
L = addbc(L, Er, 0);

%%
prefs = cheboppref;
prefs.discretization = @colloc2;
e = eigs(L, 6, prefs);
tol = 1e-10;
pass(1) = norm(e + (1:6)'.^2, inf) < tol;

%%
prefs.discretization = @ultraS;
e = eigs(L, 6, 0, prefs);
tol = 1e-10;
pass(2) = norm(e + (1:6)'.^2, inf) < tol;


%%
prefs.discretization = @colloc1;
e = eigs(L, 6, 0, prefs);
tol = 8e-10;
pass(3) = norm(e + (1:6)'.^2, inf) < tol;

end

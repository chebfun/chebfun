function pass = test_eigenvalues(pref)
% TAD, 10 Jan 2014

if ( nargin == 0 )
    pref = chebpref();
end

dom = [-pi/2, pi/2];
D2 = operatorBlock.diff(dom, 2);
E = functionalBlock.eval(dom);
El = E(dom(1));
Er = E(dom(end));
L = linop(D2);
L = addbc(L, El, 0);
L = addbc(L, Er, 0);

%%
L.prefs.discretization = @colloc2;
e = eigs(L, 6);
tol = 1e-10;
pass(1) = norm(e + (1:6)'.^2, inf) < tol;

%%
L.prefs.discretization = @ultraS;
e = eigs(L, 6, 0);
tol = 1e-10;
pass(2) = norm(e + (1:6)'.^2, inf) < tol;

end

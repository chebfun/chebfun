function pass = test_eigsRayleigh(pref)

if ( nargin < 1 )
    pref = cheboppref();
end

tol = 1e2*pref.bvpTol;

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
[v, d] = eigs(L, 6, prefs, 'rayleigh');
e = diag(d);
err(1) = norm(e - e_true, inf);
% check that we actually computed eigenfunctions
err(2) = norm(L*v-v*d);

%%
prefs.discretization = @ultraS;
[V, D] = eigs(L, 6, 0, prefs, 'rayleigh');
e = diag(D);
err(3) = norm(e - e_true, inf);
err(4) = norm(L*V-V*D);

%%
prefs.discretization = @chebcolloc1;
[V, D] = eigs(L, 6, prefs, 'rayleigh');
e = diag(D);
err(5) = norm(e - e_true, inf);
err(6) = norm(L*V-V*D);

%%
pass = err < tol;

end

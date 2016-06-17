function pass = test_vertcat(pref)

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) cos(x.*y.*z));
F = [f; f];
G = [f; f; f];
H = [F; f];
K = [f; F]; 
Fc = F.components;
Gc = G.components;
pass(1) = norm(Fc{1} - f) < tol;

pass(2) = norm(Fc{2} - f) < tol;

pass(3) = norm(Gc{1} - f) < tol;

pass(4) = norm(Gc{2} - f) < tol;

pass(5) = norm(Gc{3} - f) < tol;

pass(6) = norm(G - H) < tol;

pass(7) = norm(G - K) < tol;

f = chebfun3(@(x,y,z) cos(x.*y.*z), [-3 2 -1 2 -1 1]);
F = [f; f];
G = [f; f; f];
H = [F; f];
K = [f; F];
Fc = F.components;
Gc = G.components;
pass(8) = norm(Fc{1} - f) < tol;

pass(9) = norm(Fc{2} - f) < tol;

pass(10) = norm(Gc{1} - f) < tol;

pass(11) = norm(Gc{2} - f) < tol;

pass(12) = norm(Gc{3} - f) < tol;

pass(13) = norm(G - H) < tol;

pass(14) = norm(G - K) < tol;

end
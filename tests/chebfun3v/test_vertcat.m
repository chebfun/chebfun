function pass = test_vertcat(pref)

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.cheb3Prefs.chebfun3eps;
j = 1;

f = chebfun3(@(x,y,z) cos(x.*y.*z));
F = [f; f];
G = [f; f; f];
H = [F; f];
K = [f; F]; 
Fc = F.components;
Gc = G.components;
pass(j) = norm(Fc{1} - f) < tol; 
j = j + 1;

pass(j) = norm(Fc{2} - f) < tol; 
j = j + 1; 

pass(j) = norm(Gc{1} - f) < tol; 
j = j + 1; 

pass(j) = norm(Gc{2} - f) < tol; 
j = j + 1; 

pass(j) = norm(Gc{3} - f) < tol; 
j = j + 1; 

pass(j) = norm(G - H) < tol; 
j = j + 1; 

pass(j) = norm(G - K) < tol; 
j = j + 1; 

f = chebfun3(@(x,y,z) cos(x.*y.*z), [-3 2 -1 2 -1 1]);
F = [f; f];
G = [f; f; f];
H = [F; f];
K = [f; F];
Fc = F.components;
Gc = G.components;
pass(j) = norm(Fc{1} - f) < tol;
j = j + 1; 

pass(j) = norm(Fc{2} - f) < tol; 
j = j + 1;

pass(j) = norm(Gc{1} - f) < tol; 
j = j + 1;

pass(j) = norm(Gc{2} - f) < tol; 
j = j + 1;

pass(j) = norm(Gc{3} - f) < tol; 
j = j + 1;

pass(j) = norm(G - H) < tol; 
j = j + 1;

pass(j) = norm(G - K) < tol; 
j = j + 1; 

end
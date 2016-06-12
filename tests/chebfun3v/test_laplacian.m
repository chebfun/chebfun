function pass = test_laplacian(pref)
% Test LAPLACIAN 

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 1e2*pref.cheb3Prefs.chebfun3eps;

% Check definition: 
F = chebfun3v(@(x,y,z) cos(x), @(x,y,z) sin(y), @(x,y,z) exp(z));
lapF = laplacian(F);
lapF1 = diff(F(1), 2, 1) + diff(F(1), 2, 2) + diff(F(1), 2, 3);
lapF2 = diff(F(2), 2, 1) + diff(F(2), 2, 2) + diff(F(2), 2, 3);
lapF3 = diff(F(3), 2, 1) + diff(F(3), 2, 2) + diff(F(3), 2, 3);
pass(1) = norm(lapF1 - lapF(1)) < tol;
pass(2) = norm(lapF2 - lapF(2)) < tol;
pass(3) = norm(lapF3 - lapF(3)) < tol;

% Check definition: 
F = chebfun3v(@(x,y,z) cos(x), @(x,y,z) sin(y), @(x,y,z) x.*y + y.^2);
lapF = laplacian(F);
lapF1 = diff(F(1), 2, 1) + diff(F(1), 2, 2) + diff(F(1), 2, 3);
lapF2 = diff(F(2), 2, 1) + diff(F(2), 2, 2) + diff(F(2), 2, 3);
lapF3 = diff(F(3), 2, 1) + diff(F(3), 2, 2) + diff(F(3), 2, 3);
pass(4) = norm(lapF1 - lapF(1)) < tol;
pass(5) = norm(lapF2 - lapF(2)) < tol;
pass(6) = norm(lapF3 - lapF(3)) < tol;
end
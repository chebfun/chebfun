function pass = test_laplacian(pref)
% Test LAPLACIAN 

if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 50*pref.cheb3Prefs.chebfun3eps;

% Check definition: 
F = chebfun3(@(x,y,z) cos(x.*z)+sin(y.*z));
lapF = laplacian(F);
lapF1 = diff(F, 2, 1) + diff(F, 2, 2) + diff(F, 2, 3);
pass(1) = norm(lapF1 - lapF) < tol;

% Check definition: 
F = chebfun3(@(x,y,z) cos(x) + y.*z + z.^2);
lapF = laplacian(F);
lapF1 = diff(F, 2, 1) + diff(F, 2, 2) + diff(F, 2, 3);
pass(2) = norm(lapF1 - lapF) < tol;

end
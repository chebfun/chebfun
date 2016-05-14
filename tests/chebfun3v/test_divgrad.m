function pass = test_divgrad(pref)
% Test DIVGRAD = LAP
if ( nargin == 0 )
    pref = chebfunpref; 
end
tol = 50*pref.cheb3Prefs.chebfun3eps;

% Check definition: 
F = chebfun3v(@(x,y,z) cos(x), @(x,y,z) sin(y), @(x,y,z) exp(z));
divgradF = diffx(F(1), 2) + diffy(F(2), 2) + diffz(F(3), 2);
pass(1) = norm(divgradF - divgrad(F)) < tol;

end
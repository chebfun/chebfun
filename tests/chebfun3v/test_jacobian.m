function pass = test_jacobian(pref)
% Test Jacobian

if ( nargin == 0 ) 
    pref = chebfunpref;
end
tol = 100*pref.cheb3Prefs.chebfun3eps;

% Check definition:
F = chebfun3v(@(x,y,z) cos(x), @(x,y,z) sin(y), @(x,y,z) exp(z));
Fx = diffx(F);
Fy = diffy(F);
Fz = diffz(F);
jacF = Fx(1).*Fy(2).*Fz(3);

pass(1) = norm(jacF - jacobian(F) ) < tol;

end
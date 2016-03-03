function pass = test_jacobian( pref ) 
% Test Jacobian

if ( nargin == 0 ) 
    pref = chebfunpref; 
end
tol = 100*pref.cheb2Prefs.chebfun2eps;


% Check definition:
F = chebfun2v(@(x,y) cos(x), @(x,y) sin(y));
Fx = diffx(F); 
Fy = diffy(F); 
jacF = Fx(1).*Fy(2) - Fx(2).*Fy(1); 

pass(1) = ( norm(jacF - jacobian(F) ) < tol );

end

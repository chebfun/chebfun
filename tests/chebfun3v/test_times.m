function pass = test_times(pref)
% Check the Chebfun3v constructor for simple arithmetic operations.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
F = [f; f]; 
G = [2*f; 2*f]; 
H = [f/2; f/2];
K = [f.^2; f.^2];

pass(1) = norm(2*F - G) < tol; 

pass(2) = norm(2.*F - G) < tol;

pass(3) = norm(F*2 - G) < tol;

pass(4) = norm(F.*2 - G) < tol;

pass(5) = norm(F/2 - H) < tol;

pass(6) = norm(F./2 - H) < tol;

pass(7) = norm(2.\F - H) < tol;

pass(8) = norm(2\F - H) < tol;

pass(9) = norm(F.^2 - K) < tol;

pass(10) = norm(F.*F - K) < tol;

end
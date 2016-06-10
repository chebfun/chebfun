function pass = test_times( pref ) 
% Check the Chebfun3v constructor for simple arithmetic operations.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.cheb3Prefs.chebfun3eps;
j = 1;

f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
F = [f ; f]; 
G = [2*f ; 2*f]; 
H = [f/2 ; f/2];
K = [f.^2 ; f.^2]; 

pass(j) = norm(2*F - G) < tol; 
j = j + 1; 

pass(j) = norm(2.*F - G) < tol; 
j= j + 1; 

pass(j) = norm(F*2 - G) < tol; 
j = j + 1; 

pass(j) = norm(F.*2 - G) < tol; 
j = j + 1; 

pass(j) = norm(F/2 - H) < tol; 
j = j + 1; 

pass(j) = norm(F./2 - H) < tol; 
j = j + 1; 

pass(j) = norm(2.\F - H) < tol; 
j = j + 1; 

pass(j) = norm(2\F - H) < tol; 
j = j + 1; 

pass(j) = norm(F.^2 - K) < tol; 
j = j + 1; 

pass(j) = norm(F.*F - K) < tol; 
j = j + 1; 

end
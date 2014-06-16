function pass = test_times( pref ) 
% Check the Chebfun2v constructor for simple arithmetic operations.
% Alex Townsend, March 2013.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.eps; 
j = 1;

f = chebfun2(@(x,y) cos(x.*y)); 
F = [f ; f]; 
G = [2*f ; 2*f]; 
H = [f/2 ; f/2];
K = [f.^2 ; f.^2]; 


pass(j) = norm( 2*F - G ) < tol; j = j + 1; 
pass(j) = norm( 2.*F - G ) < tol; j = j + 1; 
pass(j) = norm( F*2 - G ) < tol; j = j + 1; 
pass(j) = norm( F.*2 - G ) < tol; j = j + 1; 
pass(j) = norm( F/2 - H ) < tol; j = j + 1; 
pass(j) = norm( F./2 - H ) < tol; j = j + 1; 
pass(j) = norm( 2.\F - H ) < tol; j = j + 1; 
pass(j) = norm( 2\F - H ) < tol; j = j + 1; 
pass(j) = norm( F.^2 - K ) < tol; j = j + 1; 
pass(j) = norm( F.*F - K ) < tol; j = j + 1; 


end
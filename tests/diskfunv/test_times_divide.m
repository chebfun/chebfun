function pass = test_times_divide( pref ) 
% Check the times and divide operations.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.cheb2Prefs.chebfun2eps;
j = 1;

f = diskfun(@(x,y) sin(y) + cos(2*x.*y)); 
F = [f ; f]; 
G = [2*f ; 2*f]; 
H = [f/2 ; f/2 ];
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
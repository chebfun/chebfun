function pass = test_get( pref ) 
% Test GET.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.cheb2Prefs.chebfun2eps;
j = 1; 

f = diskfun(@(x,y) 1 + sin(pi*x.*y) + sin(pi*x));
g = diskfun(@(x,y) cos(1-x.*y)); 
F = diskfunv(f,g);
Fc = F.components;
G = F.';
pass(j) = norm( Fc{1} - f ) < tol; j = j + 1; 
pass(j) = norm( Fc{2} - g ) < tol; j = j+1; 
pass(j) = norm( 2 - F.nComponents ) < tol; j = j + 1; 
pass(j) = F.isTransposed==0; j = j + 1; 
pass(j) = G.isTransposed==1;


end
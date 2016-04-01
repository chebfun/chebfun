function pass = test_vertcat( pref ) 
% Test for checking that spherefunv objects work correctly

% Testing chebfun2v objects with three components. 
if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.cheb2Prefs.chebfun2eps;
j = 1;

f = spherefun(@(x,y,z) cos(2*pi*x.*y.*z)); 
g = spherefun(@(x,y,z) sin(2*pi*x.*y.*z)); 
h = spherefun(@(x,y,z) 1+ x.*y.*z); 
F = [f ; g; h]; 
Gc = F.components;
pass(j) = norm( Gc{1} - f ) < tol; j = j + 1; 
pass(j) = norm( Gc{2} - g ) < tol; j = j + 1; 
pass(j) = norm( Gc{3} - h ) < tol; j = j + 1; 

G = spherefunv(f,g,h);
Gc = G.components;
pass(j) = norm( Gc{1} - f ) < tol; j = j + 1; 
pass(j) = norm( Gc{2} - g ) < tol; j = j + 1; 
pass(j) = norm( Gc{3} - h ) < tol; j = j + 1; 
pass(j) = norm( G - F ) < tol;

end

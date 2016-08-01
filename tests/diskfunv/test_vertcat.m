function pass = test_vertcat( pref ) 
% Test for checking that diskfunv objects work correctly

% Testing chebfun2v objects with three components. 
if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.cheb2Prefs.chebfun2eps;
j = 1;

f = diskfun(@(x,y) cos(2*pi*x.*y)); 
g = diskfun(@(x,y) sin(2*pi*x.*y)); 

F = [f ; g]; 
Gc = F.components;
pass(j) = norm( Gc{1} - f ) < tol; j = j + 1; 
pass(j) = norm( Gc{2} - g ) < tol; j = j + 1; 

G = diskfunv(f,g);
Gc = G.components;
pass(j) = norm( Gc{1} - f ) < tol; j = j + 1; 
pass(j) = norm( Gc{2} - g ) < tol; j = j + 1; 
pass(j) = norm( G - F ) < tol;

end

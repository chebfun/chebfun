function pass = test_vertcat( pref ) 
% A chebfun2v test for checking that chebfun2v objects with two
% components is working correctly. 

% Testing chebfun2v objects with three components. 
if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e3 * pref.eps; 
j = 1;

f = chebfun2(@(x,y) cos(x.*y)); 
F = [f ; f]; 
G = [f ; f ; f]; 
H = [ F ; f];
K = [f ; F]; 
Fc = F.components;
Gc = G.components;
pass(j) = norm( Fc{1} - f ) < tol; j = j + 1; 
pass(j) = norm( Fc{2} - f ) < tol; j = j + 1; 
pass(j) = norm( Gc{1} - f ) < tol; j = j + 1; 
pass(j) = norm( Gc{2} - f ) < tol; j = j + 1; 
pass(j) = norm( Gc{3} - f ) < tol; j = j + 1; 
pass(j) = norm( G - H ) < tol; j = j + 1; 
pass(j) = norm( G - K ) < tol; j = j + 1; 


f = chebfun2(@(x,y) cos(x.*y), [-3 2 -1 2]); 
F = [f ; f]; 
G = [f ; f ; f]; 
H = [ F ; f];
K = [f ; F]; 
Fc = F.components;
Gc = G.components;
pass(j) = norm( Fc{1} - f ) < tol; j = j + 1; 
pass(j) = norm( Fc{2} - f ) < tol; j = j + 1; 
pass(j) = norm( Gc{1} - f ) < tol; j = j + 1; 
pass(j) = norm( Gc{2} - f ) < tol; j = j + 1; 
pass(j) = norm( Gc{3} - f ) < tol; j = j + 1; 
pass(j) = norm( G - H ) < tol; j = j + 1; 
pass(j) = norm( G - K ) < tol; j = j + 1; 

end

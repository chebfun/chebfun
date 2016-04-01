function pass = test_dot( pref ) 
% Test DOT
if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 50*pref.cheb2Prefs.chebfun2eps;

% Check definition: 
F = chebfun2v(@(x,y) cos(x), @(x,y) sin(y)); 
G = chebfun2v(@(x,y) x, @(x,y) y); 
dotF1 = dot(F, G);
dotF2 = F' * G;  
pass(1) = ( norm(dotF1 - dotF2) < tol );

% Check definition: 
F = chebfun2v(@(x,y) cos(x), @(x,y) sin(y), @(x,y) x.*y); 
G = chebfun2v(@(x,y) x, @(x,y) y, @(x,y) x + y ); 
dotF1 = dot(F, G);
dotF2 = F' * G;  
pass(2) = ( norm(dotF1 - dotF2) < tol );

end
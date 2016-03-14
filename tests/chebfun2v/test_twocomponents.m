function pass = test_twocomponents( pref ) 
% A chebfun2v test for checking that chebfun2v objects with two
% components is working correctly. 

% Testing chebfun2v objects with three components. 
if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e5 * pref.cheb2Prefs.chebfun2eps;
j = 1;

% Different calls to the constructor. 
f = chebfun2(@(x,y) x ); 
F1 = chebfun2v(f,f); % from chebfun2 objects. 
F2 = chebfun2v(@(x,y) x, @(x,y) x); % handles. 
pass(j) = ( norm(F1-F2) < tol ) ; j = j+1; 

% With domains. 
d = [-1 1 -1 1];
f = chebfun2(@(x,y) x,d); 
F1 = chebfun2v(f,f,d); % from chebfun2 objects. 
F2 = chebfun2v(@(x,y) x, @(x,y) x,d); % handles. 
pass(j) = ( norm(F1-F2) < tol ) ; j = j+1; 


% Plus chebfun2v objects. 
G = chebfun2v(f,2*f); 
H = F1 + G; 
pass(j) = ( norm( H(1,1) - [2 3]' )  < tol ); j = j + 1; 
H = G + 1; 
pass(j) = ( norm( H(pi/6,1) - pi/6*[1 2]' - 1 ) <tol ); j = j + 1; 
H = G + [1 2]; 
pass(j) = ( norm( H(pi/6,1) - pi/6*[1 2]' -[1 2]' ) <tol ); j = j + 1; 

% in reverse
H = 1+ G; 
pass(j) = ( norm( H(pi/6,1) - pi/6*[1 2]' -1 ) <tol ); j = j + 1; 
H = [1 2] + G; 
pass(j) = ( norm( H(pi/6,1) - pi/6*[1 2]' -[1 2]' ) <tol ); j = j + 1; 

% .* chebfun2v objects. 
G = chebfun2v(f,2*f); 
H = G.*1; 
pass(j) = ( norm( H(pi/6,1) - G(pi/6,1) ) <tol ); j = j + 1; 
H = G.*[1 2]'; 
pass(j) = ( norm( H(pi/6,1) -  pi/6*[1 4]') <tol ); j = j + 1; 

% in reverse
H = 1.*G; 
pass(j) = ( norm( H(pi/6,1) - G(pi/6,1) ) <tol ); j = j + 1; 
H = [1 2]'.*G; 
pass(j) = ( norm( H(pi/6,1) -  pi/6*[1 4]') <tol ); j = j + 1;


% * chebfun2v objects. 
G = chebfun2v(f,2*f); 
H = G*1; 
pass(j) = ( norm( H(pi/6,1) - G(pi/6,1) ) <tol ); j = j + 1; 
H = G'*[1 2]';   % inner product
pass(j) = ( norm( H - (f+4*f)) <tol ); j = j + 1; 

% in reverse
H = 1*G; 
pass(j) = ( norm( H(pi/6,1) - G(pi/6,1) ) <tol ); j = j + 1; 
try 
    H = [1 2]'*G;  % this should fail .
catch
    pass(j) = 1 ; j=j+1; 
end

% Vector calculus identities
f = chebfun2(@(x,y) sin(x.*y)); 

pass(j) = ( norm(curl( F1 + G ) - (curl(F1) + curl(G))) < 100*tol); j=j+1; 
pass(j) = ( norm(divergence( f*G ) - (dot(G,grad(f)) + f.*divergence(G))) <100*tol); j=j+1; 
pass(j) = ( norm(diff(curl(G),1,1) + diff(curl(G),1,2)) < 10*tol); j = j + 1; 
pass(j) = ( norm(divergence(grad(f)) - laplacian(f)) < 10*tol); j = j + 1; 

end
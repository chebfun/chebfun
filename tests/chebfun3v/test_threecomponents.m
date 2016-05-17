function pass = test_threecomponents( pref ) 
% Check whether Chebfun3v works properly with three components.

% Testing chebfun3v objects with three components. 
if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e5 * pref.cheb3Prefs.chebfun3eps;
j = 1;

% Different calls to the constructor. 
f = chebfun3(@(x,y,z) x); 
F1 = chebfun3v(f, f, f); % from chebfun3 objects. 
F2 = chebfun3v(@(x,y,z) x, @(x,y,z) x, @(x,y,z) x); % handles. 
pass(j) = norm(F1-F2) < tol; 
j = j+1;

% With domains. 
dom = [-1 1 -1 1 -1 1];
f = chebfun3(@(x,y,z) x, dom); 
F1 = chebfun3v(f, f, f, dom); % from chebfun3 objects.
F2 = chebfun3v(@(x,y,z) x, @(x,y,z) x, @(x,y,z) x, dom); % handles. 
pass(j) = norm(F1-F2) < tol;
j = j+1; 

% Plus chebfun3v objects:
G = chebfun3v(f, 2*f, 3*f);
H = F1 + G; 
pass(j) = norm(H(1,1,1) - [2 3 4]')  < tol; 
j = j + 1; 

H = G + 1; 
pass(j) = norm(H(pi/6, 0.5, pi/2) - pi/6*[1 2 3]' - 1) <tol; 
j = j + 1; 

H = G + [1 2 3]; 
pass(j) = norm(H(pi/6, 0.5, pi/2) - pi/6*[1 2 3]' - [1 2 3]') <tol; 
j = j + 1; 

% in reverse
H = 1 + G; 
pass(j) = norm(H(pi/6, 0.5, pi/2) - pi/6*[1 2 3]' -1) <tol; 
j = j + 1; 

H = [1 2 3] + G; 
pass(j) = norm(H(pi/6, 0.5, pi/2) - pi/6*[1 2 3]' - [1 2 3]') <tol; 
j = j + 1; 

% .* chebfun3v objects. 
G = chebfun3v(f, 2*f, 3*f); 
H = G.*1;
pass(j) = norm(H(pi/6, 0.5, pi/2) - G(pi/6, 0.5, pi/2)) <tol; 
j = j + 1;

H = G.*[1 2 3]'; 
pass(j) = norm( H(pi/6, 0.5, pi/2) -  pi/6*[1 4 9]') <tol; 
j = j + 1; 

% in reverse
H = 1.*G;
pass(j) = norm( H(pi/6, 0.5, pi/2) - G(pi/6, 0.5, pi/2) ) <tol; 
j = j + 1; 

H = [1 2 3]'.*G; 
pass(j) = norm( H(pi/6, 0.5, pi/2) -  pi/6*[1 4 9]') <tol; 
j = j + 1;

% * chebfun3v objects. 
G = chebfun3v(f, 2*f, 3*f); 
H = G*1;
pass(j) = norm(H(pi/6, 0.5, pi/2) - G(pi/6, 0.5, pi/2) ) < tol; 
j = j + 1;

H = G'*[1 2 3]';   % inner product
pass(j) = norm(H - (f+4*f+9*f)) <tol; 
j = j + 1; 

% in reverse
H = 1*G; 
pass(j) = norm(H(pi/6, 0.5, pi/2) - G(pi/6, 0.5, pi/2)) <tol; 
j = j + 1; 
try 
    H = [1 2 3]'*G;  % this should fail .
catch
    pass(j) = 1 ; 
    j=j+1; 
end

% Vector calculus identities
f = chebfun3(@(x,y,z) sin(x.*y.*z));

pass(j) = norm(curl(F1 + G) - (curl(F1) + curl(G))) < 100*tol;
j=j+1; 

pass(j) = norm(divergence(f*G) - (dot(G, grad(f)) + f.*divergence(G))) <100*tol;
j=j+1; 

pass(j) = norm(divergence(curl(G))) < 10*tol;
j = j + 1;

pass(j) = norm(divergence(grad(f)) - laplacian(f)) < 10*tol;
j = j + 1; 

pass(j) = norm(curl(curl(G)) - (grad(divergence(G)) - laplacian(G))) < 10*tol;
j = j + 1; 

end
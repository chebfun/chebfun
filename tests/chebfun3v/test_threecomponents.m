function pass = test_threecomponents(pref)
% Check whether Chebfun3v works properly with three components.

% Testing chebfun3v objects with three components. 
if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e5 * pref.cheb3Prefs.chebfun3eps;

% Different calls to the constructor. 
f = chebfun3(@(x,y,z) x);
F1 = chebfun3v(f, f, f); % from chebfun3 objects. 
F2 = chebfun3v(@(x,y,z) x, @(x,y,z) x, @(x,y,z) x); % handles. 
pass(1) = norm(F1-F2) < tol; 

% With domains. 
dom = [-1 1 -1 1 -1 1];
f = chebfun3(@(x,y,z) x, dom); 
F1 = chebfun3v(f, f, f, dom); % from chebfun3 objects.
F2 = chebfun3v(@(x,y,z) x, @(x,y,z) x, @(x,y,z) x, dom); % handles. 
pass(2) = norm(F1-F2) < tol;

% Plus chebfun3v objects:
G = chebfun3v(f, 2*f, 3*f);
H = F1 + G;
pass(3) = norm(H(1,1,1) - [2 3 4]')  < tol; 

H = G + 1; 
pass(4) = norm(H(pi/6, 0.5, pi/2) - pi/6*[1 2 3]' - 1) <tol;

H = G + [1 2 3];
pass(5) = norm(H(pi/6, 0.5, pi/2) - pi/6*[1 2 3]' - [1 2 3]') <tol;

% in reverse
H = 1 + G;
pass(6) = norm(H(pi/6, 0.5, pi/2) - pi/6*[1 2 3]' -1) <tol;

H = [1 2 3] + G;
pass(7) = norm(H(pi/6, 0.5, pi/2) - pi/6*[1 2 3]' - [1 2 3]') <tol;

% .* chebfun3v objects. 
G = chebfun3v(f, 2*f, 3*f); 
H = G.*1;
pass(8) = norm(H(pi/6, 0.5, pi/2) - G(pi/6, 0.5, pi/2)) <tol;

H = G.*[1 2 3]'; 
pass(9) = norm( H(pi/6, 0.5, pi/2) -  pi/6*[1 4 9]') <tol;

% in reverse
H = 1.*G;
pass(10) = norm( H(pi/6, 0.5, pi/2) - G(pi/6, 0.5, pi/2) ) <tol;

H = [1 2 3]'.*G; 
pass(11) = norm( H(pi/6, 0.5, pi/2) -  pi/6*[1 4 9]') <tol;

% * chebfun3v objects. 
G = chebfun3v(f, 2*f, 3*f);
H = G*1;
pass(12) = norm(H(pi/6, 0.5, pi/2) - G(pi/6, 0.5, pi/2) ) < tol;

H = G'*[1 2 3]';   % inner product
pass(13) = norm(H - (f+4*f+9*f)) <tol;

% in reverse
H = 1*G; 
pass(14) = norm(H(pi/6, 0.5, pi/2) - G(pi/6, 0.5, pi/2)) <tol;
try 
    H = [1 2 3]'*G;  % this should fail .
catch
    pass(15) = 1;
end

% Vector calculus identities
f = chebfun3(@(x,y,z) sin(x.*y.*z));

pass(16) = norm(curl(F1 + G) - (curl(F1) + curl(G))) < 100*tol;

pass(17) = norm(divergence(f*G) - (dot(G, grad(f)) + f.*divergence(G))) ...
    <100*tol;

pass(18) = norm(divergence(curl(G))) < 10*tol;

pass(19) = norm(divergence(grad(f)) - laplacian(f)) < 10*tol;

pass(20) = norm(curl(curl(G)) - (grad(divergence(G)) - laplacian(G))) ...
    < 10*tol;

end
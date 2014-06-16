function pass = test_threecomponents( pref ) 
% A chebfun2v test for checking that chebfun2v objects with three 
% components is working correctly. 

% Testing chebfun2v objects with three components. 
if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e5 * pref.eps; 
j = 1;

% Different calls to the constructor. 
f = chebfun2(@(x,y) x ); 
F1 = chebfun2v(f,f,f); % from chebfun2 objects. 
F2 = chebfun2v(@(x,y) x, @(x,y) x, @(x,y) x); % handles. 
pass(j) = ( norm(F1-F2) < tol ) ; j = j+1; 

% With domains. 
d = [-1 1 -1 1];
f = chebfun2(@(x,y) x,d); 
F1 = chebfun2v(f,f,f,d); % from chebfun2 objects. 
F2 = chebfun2v(@(x,y) x, @(x,y) x, @(x,y) x,d); % handles. 
pass(j) = ( norm(F1-F2) < tol ) ; j = j+1; 


% Plus chebfun2v objects. 
G = chebfun2v(f,2*f,3*f); 
H = F1 + G; 
pass(j) = ( norm( H(1,1) - [2 3 4]' )  < tol ); j = j + 1; 
H = G + 1; 
pass(j) = ( norm( H(pi/6,1) - pi/6*[1 2 3]' - 1 ) <tol ); j = j + 1; 
H = G + [1 2 3]; 
pass(j) = ( norm( H(pi/6,1) - pi/6*[1 2 3]' -[1 2 3]' ) <tol ); j = j + 1; 

% in reverse
H = 1+ G; 
pass(j) = ( norm( H(pi/6,1) - pi/6*[1 2 3]' -1 ) <tol ); j = j + 1; 
H = [1 2 3] + G; 
pass(j) = ( norm( H(pi/6,1) - pi/6*[1 2 3]' -[1 2 3]' ) <tol ); j = j + 1; 

% .* chebfun2v objects. 
G = chebfun2v(f,2*f,3*f); 
H = G.*1; 
pass(j) = ( norm( H(pi/6,1) - G(pi/6,1) ) <tol ); j = j + 1; 
H = G.*[1 2 3]'; 
pass(j) = ( norm( H(pi/6,1) -  pi/6*[1 4 9]') <tol ); j = j + 1; 

% in reverse
H = 1.*G; 
pass(j) = ( norm( H(pi/6,1) - G(pi/6,1) ) <tol ); j = j + 1; 
H = [1 2 3]'.*G; 
pass(j) = ( norm( H(pi/6,1) -  pi/6*[1 4 9]') <tol ); j = j + 1;


% * chebfun2v objects. 
G = chebfun2v(f,2*f,3*f); 
H = G*1; 
pass(j) = ( norm( H(pi/6,1) - G(pi/6,1) ) <tol ); j = j + 1; 
H = G'*[1 2 3]';   % inner product
pass(j) = ( norm( H - (f+4*f+9*f)) <tol ); j = j + 1; 

% in reverse
H = 1*G; 
pass(j) = ( norm( H(pi/6,1) - G(pi/6,1) ) <tol ); j = j + 1; 
try 
    H = [1 2 3]'*G;  % this should fail .
catch
    pass(j) = 1 ; j=j+1; 
end

% Vector calculus identities
f = chebfun2(@(x,y) sin(x.*y)); 

pass(j) = ( norm(curl( F1 + G ) - (curl(F1) + curl(G))) < 100*tol); j=j+1; 
pass(j) = ( norm(divergence( f*G ) - (dot(G,[grad(f);0]) + f.*divergence(G))) <100*tol); j=j+1; 
pass(j) = ( norm(divergence(curl(G))) < 10*tol); j = j + 1; 
pass(j) = ( norm(divergence(grad(f)) - laplacian(f)) < 10*tol); j = j + 1; 
pass(j) = ( norm(curl(curl(G)) - ([grad(divergence(G));0] - laplacian(G))) < 10*tol); j = j + 1; 

% Torus example:
u = chebfun2(@(u,v) u, [0 2*pi 0 2*pi]); 
v = chebfun2(@(u,v) v, [0 2*pi 0 2*pi]);
x = (3+cos(v)).*cos(u); 
y = (3+cos(v)).*sin(u); 
z = sin(v); 
r = [x;y;z]; 
n = normal(r);
f = x;
r1 = diff(r,1,1); r2 = diff(r,1,2);
g = diff(f,1,1)*r1 + diff(f,1,2)*r2;
V = cross(n,g);

% Test diff: 
f1 = diff(V(2),1,2);
f2 = diff(r,1,2);
f = f1 * f2;
fc = f.components;
fx = fc{1}; 

% should be zero on the boundary. 
pass(j) = (norm(fx(2*pi,:))<10*tol); j = j +1; 
pass(j) = (norm(fx(0,:))<10*tol); 

end
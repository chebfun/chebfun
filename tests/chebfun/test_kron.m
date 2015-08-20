function pass = test_kron(pref) 
% Test the chebfun/kron() command. 

if ( nargin < 1 ) 
    pref = chebfunpref(); 
end 

tol = 200*pref.eps; 

%% Kronecker products resulting in CHEBFUN2 objects

% Rank 1 chebfun2 
f = chebfun(@(x) x.^2, pref);
g = chebfun(@(y) sin(y), pref);
 
h1 = chebfun2(@(x,y) x.^2.*sin(y));
h2 = chebfun2(@(x,y) y.^2.*sin(x)); 

pass(1) = (norm(h1 - kron(f', g)) < tol);
pass(2) = (norm(h2 - kron(f, g')) < tol); 
pass(3) = (norm(h1 - g*f') < tol); 
pass(4) = (norm(h2 - f*g') < tol); 

% Different domain in x and y  
d = [-2 pi -pi 2];
f = chebfun(@(x) x.^2, d(1:2), pref);
g = chebfun(@(y) sin(y), d(3:4), pref);
 
h1 = chebfun2(@(x,y) x.^2.*sin(y), d);
h2 = chebfun2(@(x,y) y.^2.*sin(x), [d(3:4) d(1:2)]); 

pass(5) = ( norm(h1 - kron(f', g)) < tol); 
pass(6) = (norm(h2 - kron(f, g')) < tol);
pass(7) = ( norm(h1 - g*f') < tol);
pass(8) = (norm(h2 - f*g') < tol); 

% Quasimatrices and rank 4 chebfun2
x = chebfun('x', [-1 1], pref);
F = [1 x x.^2 x.^4]; 
G = [1 cos(x) sin(x) x.^5]; 

h1 = chebfun2(@(x,y) 1 + x.*cos(y) + x.^2.*sin(y) + x.^4.*y.^5);
h2 = chebfun2(@(x,y) 1 + y.*cos(x) + y.^2.*sin(x) + y.^4.*x.^5);

pass(9) = (norm(h1 - kron(F', G)) < tol); 
pass(10) = (norm(h2 - kron(F, G')) < tol); 
pass(11) = (norm(h1 - G*F') < tol); 
pass(12) = (norm(h2 - F*G') < tol); 

% Different domains and quasimatrices 
x = chebfun('x',[-2 1], pref); 
y = chebfun('x',[-1 1], pref); 
d = [-2 1 -1 1]; 
F = [1 x x.^2 x.^4]; 
G = [1 cos(y) sin(y) y.^5]; 
h1 = chebfun2(@(x,y) 1 + x.*cos(y) + x.^2.*sin(y) + x.^4.*y.^5,d);
h2 = chebfun2(@(x,y) 1 + y.*cos(x) + y.^2.*sin(x) + y.^4.*x.^5,[d(3:4) d(1:2)]);

pass(13) = (norm(h1 - kron(F', G)) < tol); 
pass(14) = (norm(h2 - kron(F, G')) < tol); 
pass(15) = (norm(h1 - G*F') < tol);
pass(16) = (norm(h2 - F*G') < tol);

end

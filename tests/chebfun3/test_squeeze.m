function pass = test_squeeze(pref)
% Test SQUEEZE

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

tol = 10*pref.cheb3Prefs.chebfun3eps;
j= 1;

% This does not squeeze: 
f = chebfun3(@(x,y,z) cos(x.*y.*z)); 
g = squeeze(f);  
pass(j) = norm( f - g ) < tol;
j = j+1;

% 3D ---> 1D
f = squeeze(chebfun3(@(x,y,z) cos(x))); 
g = chebfun(@(x) cos(x)); 
pass(j) = norm(f - g) < tol; 
j = j+1;

% 3D ---> 1D
f = squeeze(chebfun3(@(x,y,z) cos(y), [-1 1 -2 3 -1 1]));
g = chebfun(@(x) cos(x), [-2 3]);
pass(j) = norm(f - g) < tol; 
j = j+1;

% 3D ---> 1D
f = squeeze(chebfun3(@(x,y,z) sin(z), [-1 1 -2 3 -pi pi]) );
g = chebfun(@(x) sin(x), [-pi pi]);
pass(j) = norm(f - g) < tol;
j = j+1;

% 3D ---> 2D
f = squeeze(chebfun3(@(x,y,z) cos(y.*z))); 
g = chebfun2(@(x,y) cos(x.*y)); 
pass(j) = norm(f - g) < tol; 
j = j+1;

% 3D ---> 2D
f = squeeze(chebfun3(@(x,y,z) cos(x.*z))); 
g = chebfun2(@(x,y) cos(x.*y)); 
pass(j) = norm(f - g) < tol; 
j = j+1;

% 3D ---> 2D
f = squeeze(chebfun3(@(x,y,z) cos(x.*y))); 
g = chebfun2(@(x,y) cos(x.*y)); 
pass(j) = norm(f - g) < tol;
j = j+1;

end
function pass = test_permute( pref ) 
% Test permute

if ( nargin == 0) 
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;
j = 1; 

% Symmetric function:
f = chebfun3(@(x,y,z) cos(x.*y.*z));
pass(j) = norm(f - permute(f, [1 2 3])) < tol; 
j = j + 1; 
pass(j) = norm(f -  permute(f, [1 3 2])) < tol;
j = j + 1; 
pass(j) = norm(f -  permute(f, [2 1 3])) < tol;
j = j + 1; 
pass(j) = norm(f -  permute(f, [2 3 1])) < tol;
j = j + 1; 
pass(j) = norm(f -  permute(f, [3 1 2])) < tol;
j = j + 1; 
pass(j) = norm(f -  permute(f, [3 2 1])) < tol;
j = j + 1; 

% Nonsymmetric function:
dom = [-3 4 -1 0 -pi pi];
f = chebfun3(@(x,y,z) cos(x).*sin(y) + z, dom);
pass(j) = norm(f - permute(f, [1 2 3])) < tol; 
j = j + 1;

g132 = chebfun3(@(x,z,y) cos(x).*sin(y) + z, [-3 4 -pi pi -1 0]);
pass(j) = norm(g132 - permute(f, [1 3 2])) < tol; 
j = j + 1;

g312 = chebfun3(@(z,x,y) cos(x).*sin(y) + z, [-pi pi -3 4  -1 0]);
pass(j) = norm(g312 - permute(f, [3 1 2])) < tol; 
j = j + 1;

g213 = chebfun3(@(y,x,z) cos(x).*sin(y) + z, [-1 0 -3 4 -pi pi]);
pass(j) = norm(g213 - permute(f, [2 1 3])) < tol; 
j = j + 1;

end
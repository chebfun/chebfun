function pass = test_permute(pref)
% Test permute

if ( nargin == 0) 
    pref = chebfunpref; 
end
tol = 1000*pref.cheb3Prefs.chebfun3eps;

% Symmetric function:
f = chebfun3(@(x,y,z) cos(x.*y.*z));
pass(1) = norm(f - permute(f, [1 2 3])) < tol; 

pass(2) = norm(f -  permute(f, [1 3 2])) < tol;

pass(3) = norm(f -  permute(f, [2 1 3])) < tol;

pass(4) = norm(f -  permute(f, [2 3 1])) < tol;

pass(5) = norm(f -  permute(f, [3 1 2])) < tol;

pass(6) = norm(f -  permute(f, [3 2 1])) < tol;

% Nonsymmetric function:
dom = [-3 4 -1 0 -pi pi];
f = chebfun3(@(x,y,z) cos(x).*sin(y) + 1i*z, dom);
pass(7) = norm(f - permute(f, [1 2 3])) < tol;

g132 = chebfun3(@(x,z,y) cos(x).*sin(y) + 1i*z, [-3 4 -pi pi -1 0]);
pass(8) = norm(g132 - permute(f, [1 3 2])) < tol;

g312 = chebfun3(@(z,x,y) cos(x).*sin(y) + 1i*z, [-pi pi -3 4  -1 0]);
pass(9) = norm(g312 - permute(f, [3 1 2])) < tol;

g213 = chebfun3(@(y,x,z) cos(x).*sin(y) + 1i*z, [-1 0 -3 4 -pi pi]);
pass(10) = norm(g213 - permute(f, [2 1 3])) < tol; 

end
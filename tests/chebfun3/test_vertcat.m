function pass = test_vertcat(pref)
% Test vertical concatenation of CHEBFUN3 objects. 

if ( nargin < 1 ) 
    pref = chebfunpref(); 
end
tol = 10*pref.cheb3Prefs.chebfun3eps;

f = chebfun3(@(x,y,z) x);
g = chebfun3(@(x,y,z) y);
h = chebfun3(@(x,y,z) z);
F = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z);
pass(1) = norm(vertcat(f, g, h) - F) < tol;

dom = [-1 0 -2 1 2 4];
f = chebfun3(@(x,y,z) x, dom); 
g = chebfun3(@(x,y,z) y, dom);
h = chebfun3(@(x,y,z) z, dom);
F = chebfun3v(@(x,y,z) x, @(x,y,z) y, @(x,y,z) z, dom);
pass(2) = norm(vertcat(f, g, h) - F) < tol; 
pass(3) = norm(vertcat(f, g, h) - [f; g; h]) < tol;
pass(4) = norm(vertcat(f) - f) < tol;

end
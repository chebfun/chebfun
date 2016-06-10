function pass = test_equiFlag(pref)
% Test 'equi' flag in Chebfun3 constructor 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 100*pref.cheb3Prefs.chebfun3eps;

% Default domain: 
dom = [-1 1 -1 1 -1 1];
ff = @(x,y,z) cos(x+y+z);
x = linspace(dom(1), dom(2), 100);
y = linspace(dom(3), dom(4), 100); 
z = linspace(dom(5), dom(6), 100);
[xx, yy, zz] = ndgrid(x, y, z);
T = ff(xx, yy, zz);
f = chebfun3(ff, dom); 
g = chebfun3(T, dom, 'equi');
pass(1) = norm(f - g) < tol ;

% A different domain: 
dom = [-1 2 -2 1 -3 0];
x = linspace(dom(1), dom(2), 100);
y = linspace(dom(3), dom(4), 100); 
z = linspace(dom(5), dom(6), 100);
[xx, yy, zz] = ndgrid(x, y, z);
T = ff(xx, yy, zz);
f = chebfun3(ff, dom); 
g = chebfun3(T , dom, 'equi' );
pass(2) = norm(f - g) < tol ; 

% Nonsymmetric function: 
dom = [-1 2 -2 1 -3 0];
ff = @(x,y,z) cos(x+2*y+3*z);
x = linspace(dom(1),dom(2),100); 
y = linspace(dom(3),dom(4),100); 
z = linspace(dom(5), dom(6), 100);
[xx, yy, zz] = ndgrid(x, y, z);
T = ff(xx, yy, zz);
f = chebfun3(ff, dom); 
g = chebfun3(T , dom, 'equi');
pass(3) = norm(f - g) < tol ; 

% Small domain; 
h = 1e-2;
dom = [1-h 1+h 1-2*h 1+2*h 1-3*h 1+3*h];
ff = @(x,y,z) cos(x+2*y+3*z);
x = linspace(dom(1),dom(2),100); 
y = linspace(dom(3),dom(4),100); 
z = linspace(dom(5), dom(6), 100);
[xx, yy, zz] = ndgrid(x, y, z);
T = ff(xx, yy, zz);
f = chebfun3(ff, dom); 
g = chebfun3(T , dom, 'equi');
pass(4) = norm(f - g) < tol ; 

end 
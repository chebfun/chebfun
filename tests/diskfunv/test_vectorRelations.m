function pass = test_vectorRelations( ) 
% Test Vector calculus relations, div, grad, curl.

tol = 3e7*chebfunpref().cheb2Prefs.chebfun2eps;

f = diskfun(@(x,y) cos((y+.1).*x));

% Div of gradient field is laplacian: 
pass(1) = norm(div(grad(f)) - laplacian(f), inf) < tol; 


% Div of curl is zero: 
pass(2) = norm(div(curl(f)), inf) < tol;

% curl of a gradient field is zero: 
pass(3) = norm(curl(grad(f)), inf) < tol; 

% test divgrad (laplacian of diskfunv)
f = diskfun(@(x,y) cos(y+x).*sin(pi*x)); 
g = diskfun(@(x,y) sin(pi*x.*y)); 
F = diskfunv(f, g); 
pass(4) = norm(divgrad(F) - (diffx(f, 2)+diffy(g, 2))) < tol; 


end
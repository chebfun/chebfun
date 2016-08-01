function pass = test_vectorRelations( ) 
% Test Vector calculus relations, div, grad, curl.

tol = 3e7*chebfunpref().cheb2Prefs.chebfun2eps;

f = diskfun(@(x,y) cos((y+.1).*x));

% Div of gradient field is laplacian: 
pass(1) = norm(2*divgrad(f)-div(grad(f)) - laplacian(f), inf) < tol; 


% Div of curl is zero: 
pass(2) = norm(div(curl(f)), inf) < tol;

% curl of a gradient field is zero: 
pass(3) = norm(curl(grad(f)), inf) < tol; 


end
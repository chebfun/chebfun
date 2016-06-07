function pass = test_vectorRelations( ) 
% Test Vector calculus relations, div, cross, curl.

tol = 3e3*chebfunpref().cheb2Prefs.chebfun2eps;

f = spherefun(@(x,y,z) cos((x+.1).*y.*z));

% Div of gradient field is laplacian: 
pass(1) = norm(div(grad(f)) - laplacian(f), inf) < tol; 

% Div of curl is zero: 
pass(2) = norm(div(curl(f)), inf) < tol;

% Vort of a gradient field is zero: 
pass(3) = norm(vort(grad(f)), inf) < tol; 

end
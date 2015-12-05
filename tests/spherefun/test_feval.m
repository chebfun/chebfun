function pass = test_feval( ) 
% Test spherefun feval. 

tol = 1e3*chebfunpref().techPrefs.eps;

f = @(x,y,z) sin(x+ y.*z) + 1;
f = redefine_function_handle( f );
g = spherefun( f );
lambda = rand; theta = rand; 
pass(1) = abs( feval(g, theta, lambda) - f(theta, lambda) ) < tol; 

% feval at vectors: 
lambda = rand(10,1); theta = rand(10,1); 
pass(2) = norm( feval(g, theta, lambda) - f(theta, lambda) ) < tol; 

% feval at vectors: 
% This breaks in feval@separableApprox.
% lambda = rand(1,10); theta = rand(1,10); 
% pass(12) = norm( feval(g, theta, lambda) - f(theta, lambda) ) < tol;
pass(3) = true;

% feval at vectors: 
lambda = rand(2,10); theta = rand(2,10); 
pass(4) = norm( feval(g, theta, lambda) - f(theta, lambda) ) < tol; 

% feval at meshgrid: 
[lambda, theta] = meshgrid( rand(3,1) ); 
pass(5) = norm( feval(g, theta, lambda) - f(theta, lambda) ) < tol; 

end 
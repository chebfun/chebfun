function pass = test_feval( ) 
% Test spherefun feval. 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

% See the random number generator.
rng(9);

f = @(x,y,z) sin(x+ y.*z) + 1;
f = redefine_function_handle(f);
g = spherefun(f);
lambda = rand; 
theta = rand; 
pass(1) = abs( g(theta, lambda) - f(theta, lambda) ) < tol;

% feval at vectors: 
lambda = rand(10, 1);
theta = rand(10, 1); 
pass(2) = norm( g(theta, lambda) - f(theta, lambda), inf ) < tol; 

% feval at row vectors: 
lambda = rand(1, 10);
theta = rand(1, 10);
pass(3) = norm(g(theta, lambda) - f(theta, lambda) ) < tol;

% feval at vectors:
lambda = rand(2, 10); 
theta = rand(2, 10);
pass(4) = norm(g(theta, lambda) - f(theta, lambda), inf ) < tol;

% feval at meshgrid: 
[lambda, theta] = meshgrid(rand(3, 1));
pass(5) = norm(g(theta, lambda) - f(theta, lambda), inf) < tol;

% feval using Cartesian coordinates
f = @(x,y,z) sin(x+ y.*z) + 1;
x = rand(1, 3); 
x = x/sqrt(sum(x.^2)); 
y = x(1, 2);
z = x(1, 3);
x = x(1, 1);
pass(6) = norm(g(x, y, z) - f(x, y, z), inf) < tol;

% feval using Cartesian coordinates at a vector of points
x = rand(10, 3); 
nrm = sqrt(sum(x.^2,2));
y = x(:, 2)./nrm; 
z = x(:, 3)./nrm; 
x = x(:, 1)./nrm;
pass(7) = norm(g(x, y, z) - f(x, y, z), inf) < tol;

% feval using Cartesian coordinates at a vector of points
x = rand(10, 3); 
y = rand(10, 3); 
z = rand(10, 3);
nrm = sqrt(x.^2 + y.^2 + z.^2); 
y = y./nrm; 
z = z./nrm; 
x = x./nrm;
pass(8) = norm(g(x, y, z) - f(x, y, z), inf) < tol;

% Check for the appropriate error flag.
try
    y = g(1,1,1);
    pass(9) = false;
catch ME
    pass(9) = strcmp(ME.identifier,'CHEBFUN:SPHEREFUN:FEVAL:pointsNotOnSphere');
end
 
end 

function f = redefine_function_handle(f)
% nargin(f) = 2, then we are already on the sphere, if nargin(f) = 3,
% then do change of variables:

if ( nargin(f) == 3 )
    % Wrap f so it can be evaluated in spherical coordinates
    f = @(lam, th) spherefun.sphf2cartf(f, lam, th, 0);
end

end
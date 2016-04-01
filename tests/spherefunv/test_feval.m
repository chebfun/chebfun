function pass = test_feval( ) 
% Test spherefunv feval. 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

% See the random number generator.
rng(9);

% feval in Cartesian coordinates at one point
f = @(x,y,z) sin(x+ y.*z) + 1;
g = @(x,y,z) cos(x+ y.*z) + 1;
h = @(x,y,z) z.*(x+ y.*z) + 1;
u = spherefunv(f,g,h);
uexact = @(x,y,z) [f(x,y,z) g(x,y,z) h(x,y,z)];
x = rand(1, 3); x = x/sqrt(sum(x.^2)); 
y = x(1, 2); z = x(1, 3); x = x(1, 1);
pass(1) = norm( u(x,y,z) - uexact(x,y,z)',inf ) < tol;

% feval in Cartesian coordinates at a vector of points
x = rand(10, 3); nrm = sqrt(sum(x.^2,2));
y = x(:, 2)./nrm; z = x(:, 3)./nrm; x = x(:, 1)./nrm;
pass(2) = norm(u(x,y,z) - uexact(x,y,z)', inf) < tol;

% feval in spherical coordinates at one point
f = @(x,y,z) sin(x+ y.*z) + 1;
g = @(x,y,z) cos(x+ y.*z) + 1;
h = @(x,y,z) z.*(x+ y.*z) + 1;
u = spherefunv(f,g,h);
f = redefine_function_handle(f);
g = redefine_function_handle(g);
h = redefine_function_handle(h);
uexact = @(lam,th) [f(lam,th) g(lam,th) h(lam,th)];
lambda = rand; 
theta = rand; 
pass(3) = norm( u(theta, lambda) - uexact(theta, lambda)',inf ) < tol;

% feval in spherical coordinates at a vector of points
lambda = rand(10, 1);
theta = rand(10, 1); 
pass(4) = norm( u(theta, lambda) - uexact(theta, lambda)',inf ) < tol;

% Check for the appropriate error flag.
try
    y = u(1,1,1);
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier,'CHEBFUN:SPHEREFUN:FEVAL:pointsNotOnSphere');
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
function pass = test_feval( ) 
% Test diskfunv feval. 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

% See the random number generator.
rng(9);

% feval in Cartesian coordinates at one point
f = @(x,y) sin(x+ y) + 1;
g = @(x,y) cos(x+ y) + 1;
u = diskfunv(f,g);
uexact = @(x,y) [f(x,y) g(x,y) ];
x = rand(1, 2); x = x/sqrt(sum(x.^2)); 
y = x(1, 2); ; x = x(1, 1);
pass(1) = norm( u(x,y) - uexact(x,y)',inf ) < tol;

% feval in Cartesian coordinates at a vector of points
x = rand(10, 2); nrm = sqrt(sum(x.^2,2));
y = x(:, 2)./nrm;  x = x(:, 1)./nrm;
pass(2) = norm(u(x,y) - uexact(x,y)', inf) < tol;

% feval in polar coordinates at one point
f = @(x,y) sin(x+ y) + 1;
g = @(x,y) cos(x+ y) + 1;
u = diskfunv(f,g);
f = redefine_function_handle(f);
g = redefine_function_handle(g);
uexact = @(t,r) [f(t,r) g(t,r)];
th = 2*rand-1; 
rad= pi*(2*rand-1); 
pass(3) = norm( u(th, rad, 'polar') - uexact(th, rad)',inf ) < tol;

% feval in polar coordinates at a vector of points
th = pi*(2*rand(10, 1)-1);
rad = 2*rand(10, 1)-1;  
pass(4) = norm( u(th, rad, 'polar') - uexact(th, rad)',inf ) < tol;


% Check for the appropriate error flag.
try
    y = u(1,2);
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.identifier,'CHEBFUN:DISKFUN:FEVAL:pointsNotOnDisk');
end
 
end 

function f = redefine_function_handle(f)
    % Wrap Cartesian f so it can be evaluated in polar coordinates
    
    f = @(th, r) diskfun.pol2cartf(f,th, r);
  
end
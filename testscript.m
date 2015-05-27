function pass = testscript( )

tol = 1e4*eps;
% Construction: 
% A few tests for construction: 

f = @(x,y,z) x.^2 + y.^2 + z.^2;
f = redefine_function_handle( f );
g = spherefun( f );
pass(1) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) exp(-cos(pi*(x+y+z)));
f = redefine_function_handle( f );
g = spherefun( f );
pass(2) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) 1-exp(x);
g = spherefun( f );
f = redefine_function_handle( f );
pass(3) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) exp(y);
g = spherefun( f );
f = redefine_function_handle( f );
pass(4) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) exp(z);
g = spherefun( f );
f = redefine_function_handle( f );
pass(5) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) cos(x.*y);
g = spherefun( f );
f = redefine_function_handle( f );
pass(6) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) sin(x.*y.*z);
g = spherefun( f );
f = redefine_function_handle( f );
pass(7) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) sin(x+ y.*z);
f = redefine_function_handle( f );
g = spherefun( f );
pass(8) = ( SampleError( f, g ) < tol ); 

f = @(x,y,z) sin(x+ y.*z) + 1;
f = redefine_function_handle( f );
g = spherefun( f );
pass(9) = ( SampleError( f, g ) < tol ); 

% feval: 
lambda = rand; theta = rand; 
pass(10) = abs( feval(g, theta, lambda) - f(theta, lambda) ) < tol; 

% feval at vectors: 
lambda = rand(10,1); theta = rand(10,1); 
pass(11) = norm( feval(g, theta, lambda) - f(theta, lambda) ) < tol; 

% feval at vectors: 
lambda = rand(1,10); theta = rand(1,10); 
pass(12) = norm( feval(g, theta, lambda) - f(theta, lambda) ) < tol; 

% feval at vectors: 
lambda = rand(2,10); theta = rand(2,10); 
pass(13) = norm( feval(g, theta, lambda) - f(theta, lambda) ) < tol; 

% feval at meshgrid: 
[lambda, theta] = meshgrid( rand(3,1) ); 
pass(14) = norm( fevalm(g, theta, lambda) - f(theta, lambda) ) < tol; 

% Compression tests
f = @(x,y,z) x + z;  % Rank 2 function
f = redefine_function_handle( f );
g = spherefun( f );
h = compress( g );
pass(15) = ( SampleError( f, h ) < tol ); 
pass(16) = size(h.BlockDiag,1) == 2;

f = @(x,y,z) exp(x);
f = redefine_function_handle( f );
g = spherefun( f );
h = compress( g );
pass(17) = ( SampleError( f, h ) < tol ); 

% Sum2 tests: 
f = @(x,y,z) x.^2 + y.^2 + z.^2; 
g = spherefun( f ); 
sum2( g ) 

end

function sample_error = SampleError( h, g ) 
m = 6; n = m;  
[x, y] = getPoints( m, n ); 
[L2, T2] = meshgrid(x, y);
F = h(L2, T2);
% hn = pi/length(g.cols);  % Have to adjust for the shift in y points. Ugly!
% approx = feval(g.cols,(y-hn)/pi-.5) * g.blockDiag * feval(g.rows,x/pi)';
% sample_error = norm( F - approx , inf );
% [C,D,R] = cdr(g);
% approx = feval(g.cols,y/pi) * g.blockDiag * feval(g.rows,x/pi)';
approx = fevalm(g,x,y);
sample_error = norm( F - approx , inf );
end

function [x, y] = getPoints( m, n )

% x = linspace(-pi, pi, 2*n+1)';  x( end ) = [ ];
% GBW: You can't just remove the pole and keep everything else equally
% spaced between -pi/2 and 3*pi/2.  The issue is that you can't keep
% both point -pi/ and 3*pi/2.
% y = linspace(-pi/2, 3*pi/2, 2*m+1)'; y( m+1 ) = [ ];

x = trigpts(2*n,[-pi pi]);
% GBW: I believe we have to sample at equally spaced points shifted by h/2
% to not sample the poles and keep an even total number of points.
y = trigpts(2*m,[-pi/2 3*pi/2]);
% y = y+0.5*pi/m; % Shift y by h/2 to avoid the poles
end

function f = redefine_function_handle( f )
% nargin( f ) = 2, then we are already on the sphere, if nargin( f ) = 3,
% then do change of variables:

if ( nargin( f ) == 3 )
    % Wrap f so it can be evaluated in spherical coordinates
    f = @(lam, th) sphf2cartf(f,lam,th);
    % Double g up.
    f = @(lam, th) sph2torus(f,lam,th);
end

end


function fdf = sph2torus(f,lam,th)

fdf = real(f(lam,th));

id = th-pi/2 > 100*eps;

if ~isempty(id) && any(id(:))
    fdf(id) = f(lam(id)-pi,pi-th(id));
end

end

function fdf = sphf2cartf(f,lam,th)

x = cos(lam).*cos(th);
y = sin(lam).*cos(th);
z = sin(th);

fdf = f(x,y,z);

end



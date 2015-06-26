function pass = testspherefun( )
% Main testing file, for now.  

pass(1) = all(test_constructor( )); 
pass(2) = all(test_feval( )); 
pass(3) = all(test_sum2( )); 
pass(4) = all(test_plus( ));
pass(5) = all(test_times( )); 
pass(6) = all(test_power( ));
pass(7) = all(test_abs( ));

end 


function pass = test_constructor( ) 
% Test the spherefun constructor 

% Get tolerance: 
tol = 1e3*chebfunpref().techPrefs.eps;

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

end

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


function pass = test_sum2( ) 
% Test spherefun sum2() command. 

tol = 1e3*chebfunpref().techPrefs.eps;

f = @(x,y,z) 1 + x + y + z; 
g = spherefun( f ); 
exact_int = 4*pi;
pass(1) = abs( sum2( g ) - exact_int ) < tol;

end 

function pass = test_plus( ) 
% Test spherefun plus() command 

tol = 1e3*chebfunpref().techPrefs.eps;

f1 = @(x,y,z) sin(pi*x.*y);  % Strictly even/pi-periodic
f2 = @(x,y,z) sin(pi*x.*z);  % Strictly odd/anti-periodic
g1 = spherefun(f1);
g2 = spherefun(f2);
gplus = g1 + g2;
fplus = redefine_function_handle( @(x,y,z) f1(x,y,z) + f2(x,y,z) );
lambda = rand; theta = rand; 
pass(1) = abs( feval(gplus, theta, lambda) - fplus(theta, lambda) ) < tol; 

f1 = @(lam,th) exp(cos(lam-1).*sin(th).*cos(th));  % Mixed symmetric terms
f2 = @(lam,th) exp(sin(lam-0.35).*sin(th).*cos(th));  % Mixed symmetric terms
g1 = spherefun(f1);
g2 = spherefun(f2);
gplus = g1 + g2;
fplus = @(lam,th) f1(lam,th) + f2(lam,th) ;
lambda = rand; theta = rand; 
pass(2) = abs( feval(gplus, theta, lambda) - fplus(theta, lambda) ) < tol; 

% Check that compression is working: 
f = spherefun(@(x,y,z) x.^2 + y.^2 + z.^2); 
r = rank( f ); 
g = f; 
for k = 1:10
    g = g + f; 
end 
pass(3) = norm( g - 11*f ) < vscale(g)*tol; 
pass(4) = ( rank( g ) - r ) == 0; 

% Check what happens with cancellation errors: 
f = spherefun(@(x,y,z) sin(x.*y.*z)); 
g = 2*f; 
pass(5) = ( norm( g - f - f ) < tol ); 

end

function pass = test_times( ) 
% Test times in SPHEREFUN 

tol = 1e3*chebfunpref().techPrefs.eps;

f = spherefun(@(x,y,z) sin(x.*y.*z)); 
pass(1) = norm( f.*f - f.^2 ) < tol; 

end 

function pass = test_power( ) 
% Test power in SPHEREFUN 

tol = 1e3*chebfunpref().techPrefs.eps;

f = spherefun(@(x,y,z) z );
g = spherefun(@(x,y,z) z.^2 );
pass(1) = norm( f.^2 - g ) < tol; 

end 

function pass = test_abs( ) 
% Test abs in SPHEREFUN 

tol = 1e3*chebfunpref().techPrefs.eps;

f = spherefun(@(x,y,z) -(x.^2 + y.^2 + z.^2) );
pass(1) = norm( abs(f) + f ) < tol; 

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
% y = trigpts(2*m,[-pi/2 3*pi/2]);
% y = trigpts(m,[0 pi]);
y = linspace(0,pi,m).';

% y = y+0.5*pi/m; % Shift y by h/2 to avoid the poles
end

function f = redefine_function_handle( f )
% nargin( f ) = 2, then we are already on the sphere, if nargin( f ) = 3,
% then do change of variables:

if ( nargin( f ) == 3 )
    % Wrap f so it can be evaluated in spherical coordinates
    f = @(lam, th) spherefun.sphf2cartf(f,lam,th,0);
%     % Double g up.
%     f = @(lam, th) sph2torus(f,lam,th);
end

end


% function fdf = sph2torus(f,lam,th)
% 
% fdf = real(f(lam,th));
% 
% id = th-pi/2 > 100*eps;
% 
% if ~isempty(id) && any(id(:))
%     fdf(id) = f(lam(id)-pi,pi-th(id));
% end
% 
% end
% 
% function fdf = sphf2cartf(f,lam,th)
% 
% x = cos(lam).*cos(th);
% y = sin(lam).*cos(th);
% z = sin(th);
% 
% fdf = f(x,y,z);
% 
% end



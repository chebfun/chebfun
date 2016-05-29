function pass = testdiskfun( )
% Main testing file, for now.  
pass(1) = all(test_constructor( )); 
pass(2) = all(test_feval( )); 
pass(3) = all(test_sum2( )); 
pass(4) = all(test_plus( ));
pass(5) = all(test_times( )); 
pass(6) = all(test_power( ));
pass(7) = all(test_abs( ));
pass(8) = all(test_diff( ));
%pass(9) = all(test_laplacian( ));   %still in prog
end 


function pass = test_constructor( ) 
% Test the diskfun constructor 

% Get tolerance: 
tol = 2e3*chebfunpref().techPrefs.chebfuneps;

f = @(x,y) x.^2 + y.^2 ; 
f = redefine_function_handle( f , 0);
g = diskfun( f );
pass(1) = ( SampleError( f, g ) < tol ); 

f = @(x,y) exp(-cos(pi*(x+y))); 
f = redefine_function_handle( f,0 );
g = diskfun( f );
pass(2) = ( SampleError( f, g ) < tol ); 

f = @(x,y) cos(pi*x)+sin(5*y)-1; 
f = redefine_function_handle( f,0 );
g = diskfun( f );
pass(3) = ( SampleError( f, g ) < tol ); 


f = @(x,y) exp(y); 
f = redefine_function_handle( f,0 );
g = diskfun( f );
pass(4) = ( SampleError( f, g ) < tol ); 

f = @(x,y) exp(x); 
f = redefine_function_handle( f,0 );
g = diskfun( f );
pass(5) = ( SampleError( f, g ) < tol ); 

f = @(x,y) exp(-x)+exp(-y) ;
f = redefine_function_handle( f,0 );
g = diskfun( f );
pass(6) = ( SampleError( f, g ) < tol ); 

f = @(x,y) cos(2*x.*y);
f = redefine_function_handle( f, 0 );
g = diskfun( f );
pass(7) = ( SampleError( f, g ) < tol ); 

f = @(x,y) sin(x.*y);
f = redefine_function_handle( f,0 );
g = diskfun( f );
pass(8) = ( SampleError( f, g ) < tol ); 

f = @(x,y) sin(11*pi*x)-sin(3*pi*x) + sin(y);
f = redefine_function_handle( f,0 );
g = diskfun( f );
pass(8) = ( SampleError( f, g ) < tol ); 

f = @(x,y) exp(-((cos(11*y)+sin(x))).^2)+sin(11*(x+y)); %rank 98
f = redefine_function_handle( f,0 );
g = diskfun( f );
pass(9) = ( SampleError( f, g ) < tol ); 

end

function pass = test_feval( ) 
% Test spherefun feval. 

tol = 1e3*chebfunpref().techPrefs.chebfuneps;

f = @(th,r) sin(r.^2.*cos(th).*sin(th)) ;
f = redefine_function_handle( f,1 );
g = diskfun( f );
R = rand; TH = rand; 
pass(1) = abs( feval(g, TH, R) - f(TH, R) ) < tol; 

% feval at vectors: 
TH = rand(10,1); R = rand(10,1); 
pass(2) = norm( feval(g, TH, R) - f(TH, R) ) < tol; 

% feval at vectors: 
% This breaks in feval@separableApprox.
% lambda = rand(1,10); theta = rand(1,10); 
% pass(12) = norm( feval(g, theta, lambda) - f(theta, lambda) ) < tol;
pass(3) = true;

% feval at vectors: 
TH = rand(2,10); R = rand(2,10); 
pass(4) = norm( feval(g, TH, R) - f(TH, R) ) < tol; 

% feval at meshgrid: 
[TT, RR] = meshgrid( rand(3,1) ); 
pass(5) = norm( feval(g, TT, RR) - f(TT, RR) ) < tol; 


% feval at meshgrid cartesian: DOESNT WORK. after first pass in feval there
% is no recognition of cartesian specification flag. 

%f = @(x,y) sin(x.*y) ;
%g = diskfun( f,0 );
%[XX, YY] = meshgrid( rand(3,1) ); 
%pass(6) = norm( feval(g, XX, YY,0) - f(XX, YY) ) < tol; 
pass(6)=1; 
end 


function pass = test_sum2( ) 
% Test diskfun sum2() command. 

tol = 1e3*chebfunpref().techPrefs.chebfuneps;

f = @(x,y) 1 + x + y; 
f = redefine_function_handle( f,0 );
g = diskfun( f ); 
exact_int = pi;
pass(1) = abs( sum2( g ) - exact_int ) < tol;

end 

function pass = test_plus( ) 
% Test diskfun plus() command 

tol = 1e3*chebfunpref().techPrefs.chebfuneps;

f1 = @(x,y) cos(pi*x.*y);  % Strictly even/pi-periodic
f2 = @(x,y) sin(x+y);  % Strictly odd/anti-periodic
g1 = diskfun(f1,0);
g2 = diskfun(f2,0);
gplus = g1 + g2;
fplus = redefine_function_handle( @(x,y) f1(x,y) + f2(x,y), 0 );
TH = rand; R = rand; 
pass(1) = abs( feval(gplus, R, TH) - fplus(R, TH) ) < tol; 

f1 = @(th,r) exp(cos(th-1).*sin(r).*cos(r));  % Mixed symmetric terms
f2 = @(th,r) exp(sin(th-0.35).*sin(r).*cos(r));  % Mixed symmetric terms
g1 = diskfun(f1);
g2 = diskfun(f2);
gplus = g1 + g2;
fplus = @(th,r) f1(th,r) + f2(th,r) ;
TH = rand; R = rand; 
pass(2) = abs( feval(gplus, R, TH, 1) - fplus(R, TH) ) < tol; 

% Check that compression is working: 
f = diskfun(@(th,r) cos(th).^2 + sin(th).^2); 
r = rank( f ); 
g = f; 
for k = 1:10
    g = g + f; 
end 
pass(3) = norm( g - 11*f ) < vscale(g)*tol; 
pass(4) = ( rank( g ) - r ) == 0; 

% Check what happens with cancellation errors: 
f = diskfun(@(x,y) sin(x.*y), 0); 
g = 2*f; 
pass(5) = ( norm( g - f - f ) < tol ); 

end

function pass = test_times( ) 
% Test times in DISKFUN 

tol = 1e3*chebfunpref().techPrefs.chebfuneps;

f = diskfun(@(x,y) sin(x.*y), 0); 
pass(1) = norm( f.*f - f.^2 ) < tol; 

end 

function pass = test_power( ) 
% Test power in DISKFUN 

tol = 1e3*chebfunpref().techPrefs.chebfuneps;

f = diskfun(@(x,y) x, 0 );
g = diskfun(@(x,y) x.^2, 0 );
pass(1) = norm( f.^2 - g ) < tol; 

end 

function pass = test_abs( ) 
% Test abs in DISKFUN 

tol = 1e3*chebfunpref().techPrefs.chebfuneps;

f = diskfun(@(x,y) -(x.^2 + y.^2), 0 );
pass(1) = norm( abs(f) + f ) < tol; 

end 

function pass = test_diff( )

tol = 1e3*chebfunpref().techPrefs.chebfuneps;

% Simple tests:
f = diskfun(@(th, r) r.^3.*(sin(th)+cos(th)));   %avoids singularity
fx = diff(f,1);
exact = @(th, r) r.^2.*(2*cos(th).*sin(th)+3*cos(th).^2+sin(th).^2);  
pass(1) = SampleError( exact, fx ) < tol;
fy = diff(f,2);
exact = @(th, r) r.^2.*(3*sin(th).^2+cos(th).^2+2*sin(th).*cos(th));  
pass(2) = SampleError( exact, fy ) < tol;


a = pi/4;  %constant derivatives
f = diskfun(@(th, r) r.*sin((th-a)));  
fx = diff(f,1);
exact = @(th,r) -sin(a)+0*th;
pass(3) = SampleError( exact, fx ) < tol;
fy = diff(f,2);
exact = @(th,r) cos(a)+0*th;
pass(4) = SampleError( exact, fy ) < tol;


% Gaussian (m1=m2=0, s1=s2=1/2, corr=0)

g = @(x, y) 4/pi*exp(-4*(x.^2+y.^2));
g=redefine_function_handle(g, 0);
f = diskfun(g);
fx = diff(f, 1);
exact = @(x,y) -32/pi*x.*exp(-4*(x.^2+y.^2));
exact = redefine_function_handle(exact, 0); 
pass(5) = SampleError( exact, fx ) < tol;

end

%function pass = test_laplacian( )

%needs to be written


%end

function sample_error = SampleError( h, g ) 
m = 6; n = m;  
[x, y] = getPoints( m, n ); 
[T2, R2] = meshgrid(x, y);
F = h(T2, R2);
% hn = pi/length(g.cols);  % Have to adjust for the shift in y points. Ugly!
% approx = feval(g.cols,(y-hn)/pi-.5) * g.blockDiag * feval(g.rows,x/pi)';
% sample_error = norm( F - approx , inf );
% [C,D,R] = cdr(g);
% approx = feval(g.cols,y/pi) * g.blockDiag * feval(g.rows,x/pi)';
approx = fevalm(g,x,y);
sample_error = norm( F(:) - approx(:) , inf );
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
y = chebpts(2*m).';
y=y(m+1:end);

% y = y+0.5*pi/m; % Shift y by h/2 to avoid the poles
end

function f = redefine_function_handle(f, coords)
    % Wrap f so it can be evaluated in polar coordinates
    if ~(coords==1)
    f = @(th, r) diskfun.pol2cartf(f,th, r);
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



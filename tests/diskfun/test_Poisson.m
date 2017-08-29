function pass = test_Poisson( ) 
% Check correctness of Poisson solver on the disk: 

tol = 2e3*chebfunpref().cheb2Prefs.chebfun2eps;
 
% Simple examples

% Example 1: use laplacian to check; use nonzero bcs 
%input is a diskfun
tru = @(t,r) exp(-r.*cos(t)-r.^2.*sin(2*t));
tru = diskfun(tru, 'polar'); 
bc = @(t,r) exp(-cos(t)-sin(2*t)); 
rhs = laplacian(tru); 
u = diskfun.poisson(rhs,bc, 100); 
pass(1) = ( norm(u - tru) < 2e4*tol ); 

%input coeffs 
rhs = coeffs2(rhs); 
u2 = diskfun.poisson(rhs,bc, 100); 
pass(2) = ( norm(u2-u) < tol) ; 

% Example 2: use function handle
bc = @(th) 0*th;              
f = @(th, r) -1 + 0*th;            
u = diskfun.poisson( f, bc, 100);
exact = diskfun(@(t,r) -.25*r.^2+.25, 'polar');
pass(3) = ( norm(u - exact) < tol ); 

%Example 3: eigenfunction
lam = 5.52007811028631;
rhs = -(lam)^2*diskfun.harmonic(0,2);
tru = diskfun.harmonic(0,2); 
u = diskfun.poisson(rhs, @(t) 0*t, 100); 
pass(4) = ( norm(tru-u) < tol) ;

end
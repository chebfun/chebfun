function pass = test_helmholtz( ) 
% Check correctness of Helmholtz solver on the disk: 

tol = 2e3*chebfunpref().cheb2Prefs.chebfun2eps;
 
% Simple examples

% Example 1: should return same result as poisson when K = 0
tru = @(t,r) exp(-r.*cos(t)-r.^2.*sin(2*t));
tru = diskfun(tru, 'polar'); 
bc = @(t,r) exp(-cos(t)-sin(2*t)); 
rhs = laplacian(tru); 
u = diskfun.poisson(rhs,bc, 100);
v = diskfun.helmholtz(rhs, 0, bc, 100); 
pass(1) = ( norm(v - tru) < 2e4*tol ); 
pass(2) = ( norm(v - u) < 2e4*tol );

%example 2: try different values of k
k = [ .05 .25, 1, pi, 7]; 
for j = 1:5
utru = diskfun(@(x,y) cos(5*(x+y)-.2)+sin(3*x.*y) );
f = lap(utru) + k(j)^2*utru;
bc = utru(:,1); 
u = diskfun.helmholtz(f, k(j), bc, 257, 256); 
pass(j+2) =  (norm(utru - u) < 5e4*tol );
end

%input coeffs
K = sqrt(2); 
utru = diskfun(@(t,r) cos(r.^5.*sin(5*t)) - (r).^2, 'polar'); 
bc = utru(:,1); 
rhs = lap(utru) + K^2*utru; 
rhs = coeffs2(rhs); 
u = diskfun.helmholtz(rhs,K,bc, 100); 
pass(8) = ( norm(u-utru) < tol) ; 

% input function handle  
K = 2; 
rhs =  @(t,r) 3*cos(r.*cos(t)) + cos(r.*sin(t));  
utru  = (1/3)*diskfun(rhs, 'polar'); 
bc = utru(:,1); 
u = diskfun.helmholtz( rhs,K, bc, 100);
pass(9) = ( norm(u - utru) < tol ); 

%Example 3: eigenfunction of laplacian
lam = 5.52007811028631^2; %eigenvalue
K = sqrt(lam +1); 
rhs = diskfun.harmonic(0,2);
u = diskfun.helmholtz(rhs,K, @(t) 0*t, 100); 
pass(10) = ( norm(u-rhs) < tol) ;

end
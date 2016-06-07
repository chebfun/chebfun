function pass = test_helmholtz( pref )
% Check that the solver works for Helmholtz equation
% Alex Townsend, March 2013. 

if ( nargin < 1 ) 
    pref = chebfunpref(); 
end 
tol = pref.cheb2Prefs.chebfun2eps;

% on [-1,1]x[-1,1]: 
mu1 = 3/pi; mu2 = pi/6; lam = mu1^2 + mu2^2; 
N = chebop2(@(u) laplacian(u) + lam * u); 
xbc = @(x) cos(mu1*x) + sin(mu1*x); 
ybc = @(y) cos(mu2*y) + sin(mu2*y); 

N.lbc = @(y) xbc(-1).*ybc(y); 
N.rbc = @(y) xbc(1).*ybc(y); 
N.ubc = @(x) xbc(x).*ybc(1); 
N.dbc = @(x) xbc(x).*ybc(-1); 

u = N \ 0; 
exact = chebfun2(@(x,y) xbc(x).*ybc(y)); 

pass(1) = ( norm(u - exact) < 100*tol );

% on rectangular domain 
mu1 = 1; mu2 = 2; lam = mu1^2 + mu2^2; 
d = [-2 3 -1/10 3]; 
N = chebop2(@(u) laplacian(u) + lam * u, d); 
xbc = @(x) cos(mu1*x) + sin(mu1*x); 
ybc = @(y) cos(mu2*y) + sin(mu2*y); 

N.lbc = @(y) xbc(d(1)).*ybc(y); 
N.rbc = @(y) xbc(d(2)).*ybc(y); 
N.dbc = @(x) xbc(x).*ybc(d(3)); 
N.ubc = @(x) xbc(x).*ybc(d(4)); 

u = N \ 0; 
exact = chebfun2(@(x,y) xbc(x).*ybc(y), d); 

pass(2) = ( norm(u - exact) < 300*tol );


% High frequency on rectangular domain 
mu1 = 12; mu2 = 10; lam = mu1^2 + mu2^2; 
d = [-2 3 -1/10 3]; 
N = chebop2(@(u) laplacian(u) + lam * u, d); 
xbc = @(x) cos(mu1*x) + sin(mu1*x); 
ybc = @(y) cos(mu2*y) + sin(mu2*y); 

N.lbc = @(y) xbc(d(1)).*ybc(y); 
N.rbc = @(y) xbc(d(2)).*ybc(y); 
N.dbc = @(x) xbc(x).*ybc(d(3)); 
N.ubc = @(x) xbc(x).*ybc(d(4)); 

u = N \ 0; 
exact = chebfun2(@(x,y) xbc(x).*ybc(y), d); 

x = chebpts(100,d(1:2));
y = chebpts(100,d(3:4)); 
[xx, yy]=meshgrid(x,y); 
pass(3) = ( norm(u(xx,yy) - exact(xx,yy),inf) < 1e9*tol );


% Higher frequency on rectangular domain 
mu1 = 12; mu2 = 50; lam = mu1^2 + mu2^2; 
d = [-2 3 -1/10 3]; 
N = chebop2(@(u) laplacian(u) + lam * u, d); 
xbc = @(x) cos(mu1*x) + sin(mu1*x); 
ybc = @(y) cos(mu2*y) + sin(mu2*y); 

N.lbc = @(y) xbc(d(1)).*ybc(y); 
N.rbc = @(y) xbc(d(2)).*ybc(y); 
N.dbc = @(x) xbc(x).*ybc(d(3)); 
N.ubc = @(x) xbc(x).*ybc(d(4)); 

u = N \ 0; 
exact = chebfun2(@(x,y) xbc(x).*ybc(y), d); 

x = chebpts(100,d(1:2));
y = chebpts(100,d(3:4)); 
[xx, yy]=meshgrid(x,y); 
pass(4) = ( norm(u(xx,yy) - exact(xx,yy),inf) < 1e9*tol );

% %%
% a = .5; k = 1; d = [a 1 0 2*pi];
% f = chebfun2(@(r,t) sin(t).*besselj(0,k*r), d); 
% N = chebop2(@(r,t,u) r.^2.*diffx(u,2) + r.*diffx(u,1) + diffy(u,2) + k^2*r.^2.*u , d); 
% N.lbc = f(d(1),:); 
% N.rbc = f(d(2),:); 
% N.ubc = f(:,d(4)); 
% N.dbc = f(:,d(3)); 
% u = N \ 0; 
% plot(u)
% norm(u - f) 
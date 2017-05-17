function pass = test_feval( ) 
% Test diskfun feval. 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

% See the random number generator.
rng(9);

f = @(t,r) sin(r.^2.*cos(t).*sin(t)) ;
g = diskfun( f, 'polar' );
R = rand; TH = rand*pi; 
pass(1) = abs( feval(g, TH, R, 'polar') - f(TH, R) ) < tol; 
pass(2) = abs ( feval(g, TH, R, 'polar') -f(TH, R) ) < tol; 

% feval at vectors: 
TH = 2*pi*rand(10,1)-pi; R = 2*rand(10,1)-1; 
pass(3) = norm( feval(g, TH, R, 'polar') - f(TH, R) ) < tol;
tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

% feval at row vectors: 
TH = pi*rand(1, 10);
R = rand(1, 10);
pass(4) = norm(g(TH, R, 'polar') - f(TH, R) ) < tol;

% feval at vectors:
TH= pi*rand(2, 10); 
R = rand(2, 10);
pass(5) = norm(g(TH, R, 'polar') - f(TH, R), inf ) < tol;

% feval at meshgrid: 
[TH, R] = meshgrid(rand(3, 1));
pass(6) = norm(g(TH, R, 'polar') - f(TH, R), inf) < tol;


% feval using Cartesian coordinates
f = @(x,y) sin(x.*y); 
[X, Y] = meshgrid((rand(3,2))-.5);
pass(7) = norm(g(X,Y)-f(X,Y), inf) < tol; 
pass(8) = norm(g(X,Y) -f(X,Y), inf) < tol;

%feval Cartesian vectors and points 
X = rand(10,1)-.5;
Y = rand(10,1)-.5;
pass(9) = norm(g(X,Y)-f(X,Y), inf) < tol; 

%feval for tensors
f = @(x,y,z) sin(x.*y).*(z.^2); 
D = diskfun(@(x,y)sin(x.*y));
g = @(x,y, z) feval(D, x,y).*z.^2;
gp = @(t,r, z) feval(D, t, r, 'polar').*z.^2; 
[tt,rr,zz] = meshgrid(linspace(-pi,pi,10),linspace(0,1,10),linspace(0,1,10));
xx = rr.*cos(tt);
yy = rr.*sin(tt);

gv = g(xx,yy,zz);
fv = f(xx,yy,zz);
gpv = gp(tt,rr,zz);
pass(10) = norm(gv(:)-fv(:), inf) < tol;
pass(11) = norm(gpv(:)-fv(:),inf) < tol; 


% Check for the appropriate error flag.
try
    y = g(1,1);
    pass(12) = false;
catch ME
    pass(12) = strcmp(ME.identifier,'CHEBFUN:DISKFUN:FEVAL:pointsNotOnDisk');
end
 
end 


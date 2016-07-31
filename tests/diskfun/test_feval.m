function pass = test_feval( ) 
% Test diskfun feval. 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

% See the random number generator.
rng(9);

f = @(t,r) sin(r.^2.*cos(t).*sin(t)) ;
g = diskfun( f, 'polar' );
R = rand; TH = rand*pi; 
pass(1) = abs( feval(g, TH, R, 'polar') - f(TH, R) ) < tol; 
pass(2) = abs ( feval(g, TH, R) -f(TH, R) ) < tol; 

% feval at vectors: 
TH = 2*pi*rand(10,1)-pi; R = 2*rand(10,1)-1; 
pass(3) = norm( feval(g, TH, R, 'polar') - f(TH, R) ) < tol;
tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;
%check that coordsetting works correctly
g.coords = 'cart';
pass(4) = norm( feval(g, TH, R, 'polar') - f(TH, R) ) < tol;
g.coords = 'polar';

% feval at row vectors: 
TH = pi*rand(1, 10);
R = rand(1, 10);
pass(5) = norm(g(TH, R) - f(TH, R) ) < tol;

% feval at vectors:
TH= pi*rand(2, 10); 
R = rand(2, 10);
pass(6) = norm(g(TH, R) - f(TH, R), inf ) < tol;

% feval at meshgrid: 
[TH, R] = meshgrid(rand(3, 1));
pass(7) = norm(g(TH, R) - f(TH, R), inf) < tol;
g.coords = 'cart';
pass(8) = norm(g(TH, R, 'polar') - f(TH, R), inf) < tol;

% feval using Cartesian coordinates
f = @(x,y) sin(x.*y); 
[X, Y] = meshgrid((rand(3,2))-.5);
pass(9) = norm(g(X,Y)-f(X,Y), inf) < tol; 
pass(10) = norm(g(X,Y, 'cart') -f(X,Y), inf) < tol;
g.coords = 'polar';
pass(11) = norm(g(X,Y, 'cart')-f(X,Y), inf) < tol; 
g.coords = 'cart';

%feval Cartesian vectors and points 
X = rand(10,1)-.5;
Y = rand(10,1)-.5;
pass(12) = norm(g(X,Y)-f(X,Y), inf) < tol; 


% Check for the appropriate error flag.
try
    y = g(1,1);
    pass(13) = false;
catch ME
    pass(13) = strcmp(ME.identifier,'CHEBFUN:DISKFUN:FEVAL:pointsNotOnDisk');
end
 
end 


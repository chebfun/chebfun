function pass = test_curl( ) 
% Test curl

tol = 1e2*chebfunpref().cheb2Prefs.chebfun2eps;

% Check curl of a diskfun is a diskfunv.
u = diskfun(@(x,y) x.*y);
u = curl(u);
pass(1) = isa(u,'diskfunv');

% Check curl of the zero function is zero
f = diskfun(@(x,y) 0*x);
pass(2) = norm(curl(f)) < tol; 

%
% Check curl of a non-zero function gives the correct result
%

% Example 1
u = curl(diskfun(@(x,y) (x.^2-y.^3)));
% Exact curl
f = diskfun(@(x,y) -3*y.^2 );
g = diskfun(@(x,y) -2*x);
exact = diskfunv(f,g);
pass(3) = norm(u-exact) < tol; 

% Example 2
u = curl(diskfun(@(x,y) cos(4*x)));
% Exact curl
f = diskfun(@(x,y) 0*x);
g = diskfun(@(x,y) 4*sin(4*x));
exact = diskfunv(f,g);
pass(4) = norm(u-exact) < 10*tol; 

% Example 3
u = curl(diskfun(@(x,y) cos(4*x.^2).*sin(y)));
% Exact curl
f = diskfun(@(x,y) cos(4*x.^2).*cos(y) );
g = diskfun(@(x,y) 8*x.*sin(4*x.^2).*sin(y));
exact = diskfunv(f,g);
pass(5) = norm(u-exact) < 30*tol; 



end
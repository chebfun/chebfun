function pass = test_curl( ) 
% Test curl

tol = 1e2*chebfunpref().cheb2Prefs.chebfun2eps;

% Check curl of a spherefun is a spherefunv.
u = spherefun(@(x,y,z) x.*y.*z);
u = curl(u);
pass(1) = isa(u,'spherefunv');

% Check curl of the zero function is zero
f = spherefun(@(x,y,z) 0*x);
pass(2) = norm(curl(f)) < tol; 

%
% Check curl of a non-zero function gives the correct result
%

% Example 1
u = curl(spherefun(@(x,y,z) (x-y).*z));
% Exact curl
f = spherefun(@(x,y,z) (x-y).*y + z.^2);
g = spherefun(@(x,y,z) (y-x).*x + z.^2);
h = spherefun(@(x,y,z) -(x+y).*z);
exact = spherefunv(f,g,h);
pass(3) = norm(u-exact) < tol; 

% Example 2
u = curl(spherefun(@(x,y,z) cos(4*z)));
% Exact curl
f = spherefun(@(x,y,z) -4*y.*sin(4*z));
g = spherefun(@(x,y,z) 4*x.*sin(4*z));
h = spherefun(@(x,y,z) 0*x);
exact = spherefunv(f,g,h);
pass(4) = norm(u-exact) < tol; 

% Example 3
u = curl(spherefun(@(x,y,z) cos(4*x)));
% Exact curl
f = spherefun(@(x,y,z) 0*x);
g = spherefun(@(x,y,z) -4*z.*sin(4*x));
h = spherefun(@(x,y,z) 4*y.*sin(4*x));
exact = spherefunv(f,g,h);
pass(5) = norm(u-exact) < 100*tol; 

% Example 4
u = curl(spherefun(@(x,y,z) cos(4*y)));
% Exact curl
f = spherefun(@(x,y,z) 4*z.*sin(4*y));
g = spherefun(@(x,y,z) 0*x);
h = spherefun(@(x,y,z) -4*x.*sin(4*y));
exact = spherefunv(f,g,h);
pass(6) = norm(u-exact) < 100*tol; 

end
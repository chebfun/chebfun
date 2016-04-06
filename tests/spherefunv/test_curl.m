function pass = test_curl( ) 
% Test curl

tol = 2e3*chebfunpref().cheb2Prefs.chebfun2eps;

% Check curl of an empty spherefunv is an empty spherefunv.
u = spherefunv;
f = curl(u);
pass(1) = isempty(f) & isa(f,'spherefunv');

% Check curl of a spherefunv is a spherefunv.
u = spherefunv.unormal;
f = curl(u);
pass(2) = isa(f,'spherefunv');

% Check curl of the zero field is zero
f = spherefun(@(x,y,z) 0*x);
u = spherefunv(f,f,f);
pass(3) = norm(curl(u)) < tol; 

%
% Check curl of a non-zero field gives the correct result
%

% Example 1
f = spherefun(@(x,y,z) x);
g = spherefun(@(x,y,z) y);
h = spherefun(@(x,y,z) z);
u = spherefunv(f,g,h);
curlu = curl(u);
% Exact curl is zero
exact = spherefun(@(x,y,z) 0*x, @(x,y,z) 0*x, @(x,y,z) 0*x);
pass(4) = norm(curlu-exact) < tol; 

% Example 2
f = spherefun(@(x,y,z) cos(4*y));
g = spherefun(@(x,y,z) cos(4*z));
h = spherefun(@(x,y,z) sin(4*x));
u = spherefunv(f,g,h);
curlu = curl(u);
% Exact curl
f = spherefun(@(x,y,z) -4*x.*y.*cos(4*x)+4*(1-z.^2).*sin(4*z));
g = spherefun(@(x,y,z) 4*y.*z.*sin(4*y)-4*(1-x.^2).*cos(4*x));
h = spherefun(@(x,y,z) 4*x.*z.*sin(4*z)+4*(1-y.^2).*sin(4*y));
exact = spherefunv(f,g,h);
pass(5) = norm(curlu-exact) < 100*tol; 

end
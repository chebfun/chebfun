function pass = test_curl( ) 
% Test curl

tol = 1e2*chebfunpref().cheb2Prefs.chebfun2eps;

% Check vorticity of an empty spherefunv is an empty spherefun.
u = spherefunv;
f = curl(u);
pass(1) = isempty(f) & isa(f,'spherefunv');

% Check curl of a spherefun is a spherefunv.
u = spherefun(@(x,y,z) x.*y.*z);
u = curl(u);
pass(2) = isa(u,'spherefunv');

% Check curl of the zero field is zero
f = spherefun(@(x,y,z) 0*x);
u = spherefunv(f,f,f);
pass(3) = norm(curl(u)) < tol; 

% Check curl of the zero function is zero
f = spherefun(@(x,y,z) 0*x);
pass(4) = norm(curl(f)) < tol; 

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
pass(5) = norm(u-exact) < tol; 

% Example 2
u = curl(spherefun(@(x,y,z) cos(4*z)));
% Exact curl
f = spherefun(@(x,y,z) -4*y.*sin(4*z));
g = spherefun(@(x,y,z) 4*x.*sin(4*z));
h = spherefun(@(x,y,z) 0*x);
exact = spherefunv(f,g,h);
pass(6) = norm(u-exact) < tol; 

% Example 3
u = curl(spherefun(@(x,y,z) cos(4*x)));
% Exact curl
f = spherefun(@(x,y,z) 0*x);
g = spherefun(@(x,y,z) -4*z.*sin(4*x));
h = spherefun(@(x,y,z) 4*y.*sin(4*x));
exact = spherefunv(f,g,h);
pass(7) = norm(u-exact) < 100*tol; 

% Example 4
u = curl(spherefun(@(x,y,z) cos(4*y)));
% Exact curl
f = spherefun(@(x,y,z) 4*z.*sin(4*y));
g = spherefun(@(x,y,z) 0*x);
h = spherefun(@(x,y,z) -4*x.*sin(4*y));
exact = spherefunv(f,g,h);
pass(8) = norm(u-exact) < 100*tol; 

% 
% % Example 2
% f = spherefun(@(x,y,z) 4*x.*z.*sin(4*z));
% g = spherefun(@(x,y,z) 4*y.*z.*sin(4*z));
% h = spherefun(@(x,y,z) -4*(1-z.^2).*sin(4*z));
% u = spherefunv(f,g,h);
% divu = div(u);
% % Exact divergence
% exact = spherefun(@(x,y,z) -8*(2*(1-z.^2).*cos(4*z) - z.*sin(4*z)));
% pass(5) = norm(divu-exact, inf) < 100*tol; 
% 
% % Example 3
% f = spherefun(@(x,y,z) 4*(x.^2-1).*sin(4*x));
% g = spherefun(@(x,y,z) 4*x.*y.*sin(4*x));
% h = spherefun(@(x,y,z) 4*x.*z.*sin(4*x));
% u = spherefunv(f,g,h);
% divu = div(u);
% % Exact divergence
% exact = spherefun(@(x,y,z) -8*(2*(1-x.^2).*cos(4*x) - x.*sin(4*x)));
% pass(6) = norm(divu-exact, inf) < 100*tol; 
% 
% % Example 4
% f = spherefun(@(x,y,z) 4*x.*y.*sin(4*y));
% g = spherefun(@(x,y,z) 4*(y.^2-1).*sin(4*y));
% h = spherefun(@(x,y,z) 4*y.*z.*sin(4*y));
% u = spherefunv(f,g,h);
% divu = div(u);
% % Exact divergence
% exact = spherefun(@(x,y,z) -8*(2*(1-y.^2).*cos(4*y) - y.*sin(4*y)));
% pass(7) = norm(divu-exact, inf) < 100*tol; 

end
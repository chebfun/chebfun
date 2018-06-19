function pass = test_grad( ) 
% Test curl

tol = 1e2*chebfunpref().cheb2Prefs.chebfun2eps;

% Check gradient of an empty spherefun is an empty spherefunv.
u = spherefun;
f = grad(u);
pass(1) = isempty(f) & isa(f,'spherefunv');

% Check gradient of the zero function is zero
f = spherefun(@(x,y,z) 0*x);
pass(2) = norm(grad(f)) < tol; 

%
% Check gradient of a non-zero function gives the correct result
%

% Example 1
u = grad(spherefun(@(x,y,z) (x-y).*z));
% Exact gradient
f = spherefun(@(x,y,z) (1+2*x.*(y-x)).*z);
g = spherefun(@(x,y,z) (-1+2*y.*(y-x)).*z);
h = spherefun(@(x,y,z) -(x-y).*(2*z.^2 - 1));
exact = spherefunv(f,g,h);
pass(3) = norm(u-exact) < tol; 

% Example 2
u = grad(spherefun(@(x,y,z) cos(4*z)));
% Exact gradient
f = spherefun(@(x,y,z) 4*x.*z.*sin(4*z));
g = spherefun(@(x,y,z) 4*y.*z.*sin(4*z));
h = spherefun(@(x,y,z) -4*(1-z.^2).*sin(4*z));
exact = spherefunv(f,g,h);
pass(4) = norm(u-exact) < tol; 

% Example 3
u = grad(spherefun(@(x,y,z) cos(4*x)));
% Exact gradient
f = spherefun(@(x,y,z) 4*(x.^2-1).*sin(4*x));
g = spherefun(@(x,y,z) 4*x.*y.*sin(4*x));
h = spherefun(@(x,y,z) 4*x.*z.*sin(4*x));
exact = spherefunv(f,g,h);
pass(5) = norm(u-exact) < 20*tol; 

% Example 4
u = grad(spherefun(@(x,y,z) cos(4*y)));
% Exact gradient
f = spherefun(@(x,y,z) 4*x.*y.*sin(4*y));
g = spherefun(@(x,y,z) 4*(y.^2-1).*sin(4*y));
h = spherefun(@(x,y,z) 4*y.*z.*sin(4*y));
exact = spherefunv(f,g,h);
pass(6) = norm(u-exact) < 10*tol; 

% Check that the grad and gradient give the same result.
gradientu = gradient(spherefun(@(x,y,z) cos(4*y)));
pass(7) = isequal(u,gradientu);

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
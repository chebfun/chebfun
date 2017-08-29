function pass = test_div( ) 
% Test divergence

tol = 1e2*chebfunpref().cheb2Prefs.chebfun2eps;

% Check divergence of an empty spherefunv is an empty spherefun.
u = spherefunv;
f = div(u);
pass(1) = isempty(f) & isa(f,'spherefun');

% Check divergence of a spherefunv is a spherefun.
u = spherefunv.unormal;
f = div(u);
pass(2) = isa(f,'spherefun');

% Check divergence of the zero field is zero
f = spherefun(@(x,y,z) 0*x);
u = spherefunv(f,f,f);
pass(3) = norm(div(u), inf) < tol; 

%
% Check divergence of a non-zero field gives the correct result
%

% Example 1
f = spherefun(@(x,y,z) (1+2*x.*(y-x)).*z);
g = spherefun(@(x,y,z) (-1+2*y.*(y-x)).*z);
h = spherefun(@(x,y,z) -(x-y).*(2*z.^2 - 1));
u = spherefunv(f,g,h);
divu = div(u);
% Exact divergence
exact = spherefun(@(x,y,z) -6*(x-y).*z);
pass(4) = norm(divu-exact, inf) < tol; 

% Example 2
f = spherefun(@(x,y,z) 4*x.*z.*sin(4*z));
g = spherefun(@(x,y,z) 4*y.*z.*sin(4*z));
h = spherefun(@(x,y,z) -4*(1-z.^2).*sin(4*z));
u = spherefunv(f,g,h);
divu = div(u);
% Exact divergence
exact = spherefun(@(x,y,z) -8*(2*(1-z.^2).*cos(4*z) - z.*sin(4*z)));
pass(5) = norm(divu-exact, inf) < 100*tol; 

% Example 3
f = spherefun(@(x,y,z) 4*(x.^2-1).*sin(4*x));
g = spherefun(@(x,y,z) 4*x.*y.*sin(4*x));
h = spherefun(@(x,y,z) 4*x.*z.*sin(4*x));
u = spherefunv(f,g,h);
divu = div(u);
% Exact divergence
exact = spherefun(@(x,y,z) -8*(2*(1-x.^2).*cos(4*x) - x.*sin(4*x)));
pass(6) = norm(divu-exact, inf) < 100*tol; 

% Example 4
f = spherefun(@(x,y,z) 4*x.*y.*sin(4*y));
g = spherefun(@(x,y,z) 4*(y.^2-1).*sin(4*y));
h = spherefun(@(x,y,z) 4*y.*z.*sin(4*y));
u = spherefunv(f,g,h);
divu = div(u);
% Exact divergence
exact = spherefun(@(x,y,z) -8*(2*(1-y.^2).*cos(4*y) - y.*sin(4*y)));
pass(7) = norm(divu-exact, inf) < 200*tol; 

% Check that the div and divergence give the same result.
divergenceu = divergence(u);
pass(8) = isequal(divu,divergenceu);

end
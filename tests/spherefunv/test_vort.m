function pass = test_vort( ) 
% Test vorticity

tol = 2e3*chebfunpref().cheb2Prefs.chebfun2eps;

% Check vorticity of an empty spherefunv is an empty spherefun.
u = spherefunv;
f = vort(u);
pass(1) = isempty(f) & isa(f,'spherefun');

% Check vorticity of a spherefunv is a spherefun.
u = spherefunv.unormal;
f = vort(u);
pass(2) = isa(f,'spherefun');

% Check vorticity of the zero field is zero
f = spherefun(@(x,y,z) 0*x);
u = spherefunv(f,f,f);
pass(3) = norm(vort(u), inf) < tol; 

%
% Check vorticity of a non-zero field gives the correct result
%

% Example 1
f = spherefun(@(x,y,z) (x-y).*y + z.^2);
g = spherefun(@(x,y,z) (y-x).*x + z.^2);
h = spherefun(@(x,y,z) -(x+y).*z);
u = spherefunv(f,g,h);
vortu = vort(u);
% Exact vorticity
exact = spherefun(@(x,y,z) -6*(x-y).*z);
pass(4) = norm(vortu-exact, inf) < tol; 

% Example 2
f = spherefun(@(x,y,z) -4*y.*sin(4*z));
g = spherefun(@(x,y,z) 4*x.*sin(4*z));
h = spherefun(@(x,y,z) 0*x);
u = spherefunv(f,g,h);
vortu = vort(u);
% Exact vorticity
exact = spherefun(@(x,y,z) -8*(2*(1-z.^2).*cos(4*z) - z.*sin(4*z)));
pass(5) = norm(vortu-exact, inf) < 100*tol; 

% Example 3
f = spherefun(@(x,y,z) 0*x);
g = spherefun(@(x,y,z) -4*z.*sin(4*x));
h = spherefun(@(x,y,z) 4*y.*sin(4*x));
u = spherefunv(f,g,h);
vortu = vort(u);
% Exact vorticity
exact = spherefun(@(x,y,z) -8*(2*(1-x.^2).*cos(4*x) - x.*sin(4*x)));
pass(6) = norm(vortu-exact, inf) < 100*tol; 

% Example 4
f = spherefun(@(x,y,z) 4*z.*sin(4*y));
g = spherefun(@(x,y,z) 0*x);
h = spherefun(@(x,y,z) -4*x.*sin(4*y));
u = spherefunv(f,g,h);
vortu = vort(u);
% Exact vorticity
exact = spherefun(@(x,y,z) -8*(2*(1-y.^2).*cos(4*y) - y.*sin(4*y)));
pass(7) = norm(vortu-exact, inf) < 100*tol; 

% Check that the vort and vorticity give the same result.
vorticityu = vorticity(u);
pass(8) = isequal(vortu,vorticityu);

end
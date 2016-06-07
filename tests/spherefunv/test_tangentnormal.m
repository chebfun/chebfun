function pass = test_tangentnormal( ) 
% Test the tangent and normal functions

tol = 1e2*chebfunpref().cheb2Prefs.chebfun2eps;

% Define a vector normal to the sphere
f = spherefun(@(x,y,z) x);
g = spherefun(@(x,y,z) y);
h = spherefun(@(x,y,z) z);
u = spherefunv(f,g,h);

% Check that projetion of this vector in the normal direction just returns
% exactly this vector.
nu = normal(u);
pass(1) = norm(u - nu) < tol; 

% Check that projetion of this vector onto the tangent space of the sphere
% gives zero.
tu = tangent(u);
pass(2) = norm(tu) < tol; 

% Generate a vector field tangent the sphere
f = spherefun(@(x,y,z) cos(2*pi*x.*y.*z));
u = grad(f);

% Check that projetion of this vector onto the tangent space of the sphere
% gives this vector field back.
tu = tangent(u);
pass(3) = norm(u - tu) < tol; 

% Check that projetion of this vector in the normal direction just returns
% exactly this vector.
nu = normal(u);
pass(4) = norm(nu) < tol;

% Check that unormal returns the unit normal vector field to the sphere.
f = spherefun(@(x,y,z) x);
g = spherefun(@(x,y,z) y);
h = spherefun(@(x,y,z) z);
u = spherefunv(f,g,h);
un = spherefunv.unormal();
pass(5) = norm(u-un) < tol;

end
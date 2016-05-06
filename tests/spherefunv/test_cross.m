function pass = test_cross( ) 
% Test cross product

tol = 1e3*chebfunpref().cheb2Prefs.chebfun2eps;

% Test cross product of empty spherefunv.
f = spherefunv;
g = spherefunv;
h = cross(f,g);

pass(1) = isempty(h);

% Test that the cross product of the same vector valued function is zero.
f = spherefun(@(x,y,z) cos((x+.1).*y.*z));
u = grad(f);
h = cross(u,u);

pass(2) = norm(h) < tol;

% Test that the cross product of two vector-valued functions tangent to the
% sphere is normal to the sphere.
f = spherefun(@(x,y,z) cos((x+.1).*y.*z));
u = grad(f);
g = spherefun(@(x,y,z) sin(y.*z));
v = grad(g);
w = cross(u,v);
nrml = spherefunv(spherefun(@(x,y,z) x),spherefun(@(x,y,z) y),spherefun(@(x,y,z) z));
h = cross(nrml,w);

pass(3) = norm(h) < tol;

end
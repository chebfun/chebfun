function pass = test_dot( ) 
% Test dot product

tol = 10*chebfunpref().cheb2Prefs.chebfun2eps;

% Test dot product of empty spherefunv.
f = spherefunv;
g = spherefunv;
h = dot(f,g);
pass(1) = isempty(h);

% Test that the dot product of two orthogonal vectors is zero.
f = spherefun(@(x,y,z) cos((x+.1).*y.*z));
% Vector field tangent to the sphere
u = grad(f);
% Verctor field normal to the sphere
nrml = spherefunv(spherefun(@(x,y,z) x),spherefun(@(x,y,z) y),spherefun(@(x,y,z) z));
h = dot(u,nrml);
pass(2) = norm(h,inf) < tol;

% Test that the dot product gives the correct result 
u1 = spherefun(@(x,y,z) x.*z.*cos(2*y));
u2 = spherefun(@(x,y,z) y.*z.*sin(2*x));
u3 = spherefun(@(x,y,z) exp(x.*y.*z));
u = spherefunv(u1,u2,u3);
v1 = spherefun(@(x,y,z) x.*y);
v2 = spherefun(@(x,y,z) y.*z);
v3 = spherefun(@(x,y,z) z.*x);
v = spherefunv(v1,v2,v3);
f = dot(u,v);
% Same as the dot product
g = u1.*v1 + u2.*v2 + u3.*v3;
pass(3) = norm(g-f) < tol;

% Same as dot product 
g = u'*v;
pass(4) = norm(g-f) < tol;

end
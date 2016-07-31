function pass = test_dot( ) 
% Test dot product

tol = 10*chebfunpref().cheb2Prefs.chebfun2eps;

% Test dot product of empty diskfunv.
f = diskfunv;
g = diskfunv;
h = dot(f,g);
pass(1) = isempty(h);

% Check definition: 
F = chebfun2v(@(x,y) cos(x), @(x,y) sin(y)); 
G = chebfun2v(@(x,y) x, @(x,y) y); 
dotF1 = dot(F, G);
dotF2 = F' * G;  
pass(2) = ( norm(dotF1 - dotF2) < tol );

% Test with diskfun 
u1 = diskfun(@(x,y) x.*cos(2*y));
u2 = diskfun(@(x,y) y.*sin(2*x));
u = diskfunv(u1,u2);
v1 = diskfun(@(x,y) x.*y);
v2 = diskfun(@(x,y) y.^2);
v = diskfunv(v1,v2);
f = dot(u,v);
% Same as the dot product
g = u1.*v1 + u2.*v2;
pass(3) = norm(g-f) < tol;

% Same as dot product 
g = u'*v;
pass(4) = norm(g-f) < tol;

end
function pass = test_cross( ) 
% Test cross product

tol = 1e3*chebfunpref().cheb2Prefs.chebfun2eps;

% Test cross product of empty diskfunv.
f = diskfunv;
g = diskfunv;
h = cross(f,g);

pass(1) = isempty(h);

%check definition
F = diskfunv(@(x,y) cos(x), @(x,y) sin(y)); 
G = diskfunv(@(x,y) x, @(x,y) y); 
crossF = F(1) .* G(2) - F(2) .* G(1);
pass(2) = ( norm(cross(F,G) - crossF) < tol );

% Test that the cross product of the same vector valued function is zero.
f = diskfun(@(x,y) cos((x+.1).*y));
u = grad(f);
h = cross(u,u);
pass(3) = norm(h) < tol;

end
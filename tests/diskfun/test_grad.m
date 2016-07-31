function pass = test_grad( ) 
% Test curl

tol = 1e2*chebfunpref().cheb2Prefs.chebfun2eps;

% Check gradient of an empty diskfun is an empty diskfunv.
u = diskfun;
f = grad(u);
pass(1) = isempty(f) & isa(f,'diskfunv');

% Check gradient of the zero function is zero
f = diskfun(@(x,y) 0*x);
pass(2) = norm(grad(f)) < tol; 

%
% Check gradient of a non-zero function gives the correct result
%

% Example 1
u = grad(diskfun(@(x,y,z) x.^4.*(y-y.^2)) );
% Exact gradient
f = diskfun(@(x,y) 4*x.^3.*(y-y.^2));
g = diskfun(@(x,y) x.^4.*(1-2*y));
exact = diskfunv(f,g);
pass(3) = norm(u-exact) < tol; 

% Example 2
u = grad(diskfun(@(x,y) cos(4*x).*sin(x.*y)));
% Exact gradient
f = diskfun(@(x,y) -4*sin(4*x).*sin(x.*y)+cos(4*x).*cos(x.*y).*y);
g = diskfun(@(x,y) cos(4*x).*cos(x.*y).*x);
exact = diskfunv(f,g);
pass(4) = norm(u-exact) < 2*tol; 

% Example 3
u = grad(diskfun(@(x,y) exp(-3*(x.^2+(y+.2).^2)) ));
% Exact gradient
f = diskfun(@(x,y) -6.*x.*exp(-3*(x.^2+(y+.2).^2)));
g = diskfun(@(x,y) -6.*(y+.2).*exp(-3*(x.^2+(y+.2).^2)));
exact = diskfunv(f,g);
pass(5) = norm(u-exact) < 10*tol; 


% Check that the grad and gradient give the same result.
gradientu = gradient(diskfun(@(x,y) exp(-3*(x.^2+(y+.2).^2)) ));
pass(6) = isequal(u,gradientu);

end
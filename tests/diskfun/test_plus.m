function pass = test_plus( ) 
% Test spherefun plus() command 

tol = 1000*chebfunpref().cheb2Prefs.chebfun2eps;

f1 = @(x,y,z) sin(pi*x.*y);  % Strictly even/pi-periodic
f2 = @(x,y,z) sin(pi*x.*z);  % Strictly odd/anti-periodic
g1 = spherefun(f1);
g2 = spherefun(f2);
gplus = g1 + g2;
fplus = @(lam, th) spherefun.sphf2cartf(@(x,y,z) f1(x, y, z) + ...
    f2(x, y, z), lam, th, 0);
lambda = rand; 
theta = rand; 
pass(1) = abs(feval(gplus, theta, lambda) - fplus(theta, lambda)) < tol;

f1 = @(lam,th) exp(cos(lam-1).*sin(th).*cos(th));  % Mixed symmetric terms
f2 = @(lam,th) exp(sin(lam-0.35).*sin(th).*cos(th));  % Mixed symmetric terms
g1 = spherefun(f1);
g2 = spherefun(f2);
gplus = g1 + g2;
fplus = @(lam,th) f1(lam, th) + f2(lam, th);
lambda = rand; 
theta = rand;
pass(2) = abs(feval(gplus, theta, lambda) - fplus(theta, lambda)) < tol; 

% Check that compression is working: 
f = spherefun(@(x,y,z) x.^2 + y.^2 + z.^2);
r = rank(f);
g = f; 
for k = 1:10
    g = g + f;
end 
pass(3) = norm(g - 11*f, inf) < vscale(g)*tol;
pass(4) = (rank(g) - r ) == 0;

% Check what happens with cancellation errors: 
f = spherefun(@(x,y,z) sin(x.*y.*z)); 
g = 2*f;
pass(5) = ( norm(g - f - f, inf) < tol );

% Check for the case where f and g are non-zero at the poles bug 
% f+g is zero at the poles.
f = spherefun(@(x,y,z) z);
g = spherefun(@(x,y,z) z.^3);
h = f - g;
pass(6) = h.nonZeroPoles == 0;

% Check for the case where f non-zero at the poles but g is zero at the
% poles.
f = spherefun(@(x,y,z) z.*(1-z.^2));
g = spherefun(@(x,y,z) z);
h = f + g;
pass(7) = h.nonZeroPoles == 1;
h = g + f;
pass(8) = h.nonZeroPoles == 1;

end
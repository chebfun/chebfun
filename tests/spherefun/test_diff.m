function pass = test_diff( )

tol = 2e2*chebfunpref().cheb2Prefs.chebfun2eps;

% Simple tests:
f = spherefun(@(lam,th) cos(lam).*sin(th));  % x
fx = diff(f, 1);
exact = @(lam,th) 1-(cos(lam).*sin(th)).^2;  % 1-x^2 or y^2 + z^2
pass(1) = SampleError(exact, fx) < tol;
fy = diff(f, 2);
exact = @(lam,th) -cos(lam).*sin(th).*sin(lam).*sin(th);  % -x*y
pass(2) = SampleError(exact, fy) < tol;
fz = diff(f, 3);
exact = @(lam,th) -cos(lam).*sin(th).*cos(th);  % -x*z
pass(3) = SampleError(exact, fz) < tol;

f = spherefun(@(lam,th) cos(th));  % z
fx = diff(f, 1);
exact = @(lam,th) -cos(lam).*sin(th).*cos(th);  % -x*z
pass(4) = SampleError(exact, fx) < tol;
fy = diff(f, 2);
exact = @(lam,th) -sin(lam).*sin(th).*cos(th);  % -y*z
pass(5) = SampleError(exact, fy) < tol;
fz = diff(f, 3);
exact = @(lam,th) 1-cos(th).^2;  % 1-z.^2  or x^2 + y^2
pass(6) = SampleError(exact, fz) < tol;

a = 1/sqrt(2);
f = spherefun(@(lam,th) sin(lam-a).*sin(th));  % Shifted y
fx = diff(f, 1);
exact = @(lam,th) 0.25*(2*cos(lam).*cos(2*th).*sin(lam-a) - sin(2*lam-a) ...
    - 3*sin(a));
pass(7) = SampleError(exact, fx) < tol;
fy = diff(f, 2);
exact = @(lam,th) 0.25*(2*sin(lam).*cos(2*th).*sin(lam-a) + cos(2*lam-a) ...
    + 3*cos(a));
pass(8) = SampleError(exact, fy) < tol;
fz = diff(f, 3);
exact = @(lam,th) -0.5*sin(lam-a).*sin(2*th);
pass(9) = SampleError(exact, fz) < tol;

% Gaussian
lam0 = pi/sqrt(2); 
th0 = (sqrt(5) - 2)/2; 
sig2 = 1;
r2 = @(lam,th) 2*(1 - (sin(th)*sin(th0)).*cos(lam-lam0)-cos(th)*cos(th0));
f = @(lam,th) exp(-sig2 * r2(lam, th));
fx = diff(spherefun(f), 1);
exact = @(lam,th) -0.5*sig2*exp(-sig2*r2(lam,th)).*(2*cos(lam)...
    .*cos(th0).*sin(2*th)+(cos(2*lam-lam0)-3*cos(lam0)-...
    2*cos(lam).*cos(lam-lam0).*cos(2*th))*sin(th0));
pass(10) = SampleError(exact, fx) < tol;

end

function sample_error = SampleError(h, g)
m = 6; 
n = m;
[x, y] = getPoints(m, n);
[L2, T2] = meshgrid(x, y);
F = h(L2, T2);
approx = fevalm(g, x, y);
sample_error = norm(F(:) - approx(:), inf);
end

function [x, y] = getPoints(m, n)

x = trigpts(2*n,[-pi pi]);
y = linspace(0,pi,m).';

end
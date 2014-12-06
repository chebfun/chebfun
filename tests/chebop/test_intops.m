function pass = test_intops(pref)
% Test integral operators
% Toby Driscoll 28 May 2009
% Nick Hale 6 Jan 2012

% NOTE: Taken from V4 test chebop_intops.

if ( nargin == 0 )
    pref = cheboppref();
end
tol = 1e-10;

% From http://en.wikipedia.org/wiki/Integro-differential_equation
d = [0 5];
x = chebfun('x', d);
N = chebop(@(u) diff(u) + 2*u + 5*cumsum(u), d);
N.bc = @(x,u) u(0);
u = mldivide(N, 1, pref);
u_exact = 0.5*exp(-x).*sin(2*x);
err(1) = norm(u - u_exact);

% Fredholm
d = [0, 1];
x = chebfun(@(x) x, d);
K = @(x,y) sin(2*pi*(x-y));
A = chebop(@(x,u) u + fred(K,u), d);
u = x.*exp(x);
f = A*u;
err(2) = norm(u - mldivide(A, f, pref));

% Volterra
d = [0, pi];
x = chebfun(@(x) x, d);
K = @(x,y) x.*y;
A = chebop(@(x,u) u - volt(K,u), d);
f = x.^2.*cos(x) + (1-x).*sin(x);
u = mldivide(A, f, pref);
err(3) = norm( u - sin(x) );
err(4) = norm( A*u - f );

pass = err < tol;

end

